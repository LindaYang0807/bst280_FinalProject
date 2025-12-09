#!/usr/bin/env Rscript

# Exploratory analysis + differential expression + survival analyses
# GTEx breast (normal) vs TCGA BRCA (tumor) recount3 datasets.
# Outputs summary tables and plots to results/.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
  library(matrixStats)
  library(scales)
  library(patchwork)
  library(survival)
})

theme_set(theme_bw())

data_dir <- "data"
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

gtex <- readRDS(file.path(data_dir, "gtex_breast.rds"))
tcga <- readRDS(file.path(data_dir, "tcga_brca.rds"))

message("GTEx: ", nrow(gtex), " genes x ", ncol(gtex), " samples.")
message("TCGA: ", nrow(tcga), " genes x ", ncol(tcga), " samples.")
message("Assays: ", paste(assayNames(gtex), collapse = ", "))

to_tibble <- function(se) {
  as_tibble(colData(se), rownames = "sample_id")
}

strip_version <- function(x) sub("\\.\\d+$", "", x)

# Collapse duplicated gene IDs after removing version suffix.
collapse_by_gene <- function(mat, fun = c("sum", "mean")) {
  fun <- match.arg(fun)
  ids <- strip_version(rownames(mat))
  groups <- split(seq_along(ids), ids)
  out <- matrix(NA_real_, nrow = length(groups), ncol = ncol(mat),
                dimnames = list(names(groups), colnames(mat)))
  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    if (fun == "sum") {
      out[i, ] <- colSums(mat[idx, , drop = FALSE])
    } else {
      out[i, ] <- colMeans(mat[idx, , drop = FALSE])
    }
  }
  out
}

km_step_df <- function(data, time_col, event_col, group_col) {
  groups <- unique(na.omit(data[[group_col]]))
  res <- map_dfr(groups, function(g) {
    sub <- data %>% filter(.data[[group_col]] == g)
    if (nrow(sub) < 2) return(NULL)
    sf <- survfit(Surv(sub[[time_col]], sub[[event_col]]) ~ 1)
    if (length(sf$time) == 0) return(NULL)
    tibble(
      time = sf$time,
      surv = sf$surv,
      group = g
    )
  })
  res
}

# Metadata summaries
gtex_meta <- to_tibble(gtex) %>%
  transmute(
    sample_id,
    dataset = "GTEx",
    sex = gtex.sex,
    age_group = gtex.age,
    rna_integrity = gtex.smrin
  )

write_csv(count(gtex_meta, sex, name = "n"), file.path(res_dir, "gtex_sex_counts.csv"))
write_csv(count(gtex_meta, age_group, name = "n"), file.path(res_dir, "gtex_age_counts.csv"))

tcga_meta <- to_tibble(tcga) %>%
  transmute(
    sample_id,
    dataset = "TCGA",
    sample_type = tcga.gdc_cases.samples.sample_type,
    tumor_stage = tcga.gdc_cases.diagnoses.tumor_stage,
    vital_status = tcga.gdc_cases.diagnoses.vital_status,
    days_to_death = tcga.gdc_cases.diagnoses.days_to_death,
    days_to_last_follow_up = tcga.gdc_cases.diagnoses.days_to_last_follow_up,
    age_at_diagnosis = tcga.gdc_cases.diagnoses.age_at_diagnosis
  )

write_csv(count(tcga_meta, sample_type, name = "n"), file.path(res_dir, "tcga_sample_type_counts.csv"))
write_csv(count(tcga_meta, tumor_stage, name = "n"), file.path(res_dir, "tcga_stage_counts.csv"))
write_csv(count(tcga_meta, vital_status, name = "n"), file.path(res_dir, "tcga_vital_counts.csv"))

# Survival-ready table
survival_df <- tcga_meta %>%
  mutate(
    event = if_else(vital_status == "dead", 1, 0),
    time_days = coalesce(days_to_death, days_to_last_follow_up)
  ) %>%
  select(sample_id, event, time_days, tumor_stage, age_at_diagnosis)

write_csv(survival_df, file.path(res_dir, "tcga_survival_ready.csv"))

# Collapse assays to versionless gene IDs
gtex_counts <- collapse_by_gene(assay(gtex, "raw_counts"), fun = "sum")
tcga_counts <- collapse_by_gene(assay(tcga, "raw_counts"), fun = "sum")
gtex_logtpm <- collapse_by_gene(assay(gtex, "logtpm"), fun = "mean")
tcga_logtpm <- collapse_by_gene(assay(tcga, "logtpm"), fun = "mean")

common_genes <- intersect(rownames(gtex_counts), rownames(tcga_counts))
message("Common genes after version stripping: ", length(common_genes))

gtex_counts_common <- gtex_counts[common_genes, , drop = FALSE]
tcga_counts_common <- tcga_counts[common_genes, , drop = FALSE]
gtex_logtpm_common <- gtex_logtpm[common_genes, , drop = FALSE]
tcga_logtpm_common <- tcga_logtpm[common_genes, , drop = FALSE]

# Library size distributions
lib_df <- bind_rows(
  tibble(
    dataset = "GTEx",
    sample_id = colnames(gtex),
    library_size = colSums(assay(gtex, "raw_counts"))
  ),
  tibble(
    dataset = "TCGA",
    sample_id = colnames(tcga),
    library_size = colSums(assay(tcga, "raw_counts"))
  )
)

p_lib <- ggplot(lib_df, aes(x = library_size / 1e6, fill = dataset)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  scale_x_continuous("Library size (millions of counts)", labels = comma) +
  scale_y_continuous("Sample count") +
  guides(fill = "none")
ggsave(file.path(res_dir, "library_size_hist.png"), p_lib, width = 8, height = 4, dpi = 200)

# Median expression per sample
dens_df <- bind_rows(
  tibble(
    dataset = "GTEx",
    sample_id = colnames(gtex),
    logtpm_median = matrixStats::colMedians(assay(gtex, "logtpm"), na.rm = TRUE)
  ),
  tibble(
    dataset = "TCGA",
    sample_id = colnames(tcga),
    logtpm_median = matrixStats::colMedians(assay(tcga, "logtpm"), na.rm = TRUE)
  )
)

p_density <- ggplot(dens_df, aes(x = logtpm_median, fill = dataset)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous("Median log(TPM)") +
  scale_y_continuous("Density")
ggsave(file.path(res_dir, "logtpm_median_density.png"), p_density, width = 8, height = 4, dpi = 200)

# Gene overlap and PCA on top variable genes across datasets
combined_logtpm <- cbind(gtex_logtpm_common, tcga_logtpm_common)
top_var <- head(order(rowVars(combined_logtpm, na.rm = TRUE), decreasing = TRUE), 2000)
pca <- prcomp(t(combined_logtpm[top_var, ]), center = TRUE, scale. = TRUE)

tcga_group <- if ("tcga.gdc_cases.samples.sample_type" %in% colnames(colData(tcga))) {
  colData(tcga)$tcga.gdc_cases.samples.sample_type
} else {
  rep("TCGA", ncol(tcga))
}

pca_df <- tibble(
  sample_id = c(colnames(gtex), colnames(tcga)),
  dataset = c(rep("GTEx", ncol(gtex)), rep("TCGA", ncol(tcga))),
  group = c(rep("normal", ncol(gtex)), tcga_group),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)

p_pca12 <- ggplot(pca_df, aes(PC1, PC2, color = dataset, shape = group)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(GTEx = "#1b9e77", TCGA = "#d95f02")) +
  guides(shape = guide_legend(title = "Sample type"))

p_pca13 <- ggplot(pca_df, aes(PC1, PC3, color = dataset, shape = group)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(GTEx = "#1b9e77", TCGA = "#d95f02")) +
  guides(shape = "none")

ggsave(file.path(res_dir, "pca_combined.png"), p_pca12 + p_pca13, width = 10, height = 4, dpi = 200)

# Optional ComBat-seq batch adjustment (may be confounded here)
combat_counts <- NULL
if (requireNamespace("sva", quietly = TRUE)) {
  combat_counts <- tryCatch(
    {
      sva::ComBat_seq(as.matrix(cbind(gtex_counts_common, tcga_counts_common)),
                      batch = factor(c(rep("GTEx", ncol(gtex)), rep("TCGA", ncol(tcga)))),
                      group = factor(c(rep("normal", ncol(gtex)), rep("tumor", ncol(tcga))), levels = c("normal", "tumor")))
    },
    error = function(e) {
      message("ComBat_seq failed: ", e$message)
      NULL
    }
  )
} else {
  message("sva not installed; skipping ComBat_seq.")
}

# Differential expression: GTEx (normal) vs TCGA (tumor)
group <- factor(c(rep("normal", ncol(gtex)), rep("tumor", ncol(tcga))), levels = c("normal", "tumor"))
counts_combined <- cbind(gtex_counts_common, tcga_counts_common)

has_limma <- requireNamespace("limma", quietly = TRUE) && requireNamespace("edgeR", quietly = TRUE)

if (has_limma) {
  message("Using limma-voom for DE.")
  library(limma)
  library(edgeR)

  dge <- DGEList(counts = counts_combined, group = group)
  keep <- rowSums(cpm(dge) > 1) >= 10
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ group)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  de_tbl <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P") %>%
    as_tibble(rownames = "gene_id") %>%
    mutate(
      mean_logCPM_normal = rowMeans(v$E[, group == "normal", drop = FALSE]),
      mean_logCPM_tumor = rowMeans(v$E[, group == "tumor", drop = FALSE])
    )
  de_expr_mat <- v$E
} else {
  message("limma/edgeR not installed; using log-CPM + Welch t-test fallback (no voom).")
  lib_sizes <- colSums(counts_combined)
  cpm_mat <- t(t(counts_combined) / lib_sizes * 1e6)
  keep <- rowSums(cpm_mat > 1) >= 10
  cpm_mat <- cpm_mat[keep, , drop = FALSE]
  log_cpm <- log2(cpm_mat + 1)
  normal_idx <- group == "normal"
  tumor_idx <- group == "tumor"
  t_stats <- apply(log_cpm, 1, function(x) {
    res <- t.test(x[normal_idx], x[tumor_idx])
    c(stat = unname(res$statistic), p = res$p.value)
  })
  t_stats <- t(t_stats)
  logFC <- rowMeans(log_cpm[, tumor_idx, drop = FALSE]) - rowMeans(log_cpm[, normal_idx, drop = FALSE])
  de_tbl <- tibble(
    gene_id = rownames(log_cpm),
    logFC = logFC,
    AveExpr = rowMeans(log_cpm),
    t = t_stats[, "stat"],
    P.Value = t_stats[, "p"],
    adj.P.Val = p.adjust(t_stats[, "p"], method = "BH"),
    mean_logCPM_normal = rowMeans(log_cpm[, normal_idx, drop = FALSE]),
    mean_logCPM_tumor = rowMeans(log_cpm[, tumor_idx, drop = FALSE])
  )
  de_expr_mat <- log_cpm
}

write_csv(de_tbl, file.path(res_dir, "de_limma_voom.csv"))

volcano_df <- de_tbl %>%
  mutate(sig = adj.P.Val < 0.05 & abs(logFC) > 1)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = "log2 fold-change (tumor vs normal)", y = "-log10(FDR)") +
  guides(color = "none")

p_ma <- ggplot(de_tbl, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "grey70")) +
  labs(x = "Average log2 expression", y = "log2 fold-change (tumor vs normal)") +
  guides(color = "none")

ggsave(file.path(res_dir, "de_volcano.png"), p_volcano, width = 7, height = 5, dpi = 200)
ggsave(file.path(res_dir, "de_MA.png"), p_ma, width = 7, height = 5, dpi = 200)

# Co-expression network on top DE genes (default top 50 by adjusted p)
top_n_net <- 50
genes_for_net <- head(de_tbl$gene_id, top_n_net)
genes_for_net <- genes_for_net[genes_for_net %in% rownames(de_expr_mat)]
if (length(genes_for_net) >= 5) {
  cor_mat <- cor(t(de_expr_mat[genes_for_net, , drop = FALSE]), use = "pairwise.complete.obs")
  diag(cor_mat) <- NA
  cor_long <- as_tibble(as.table(cor_mat), .name_repair = "minimal")
  colnames(cor_long)[1:3] <- c("gene1", "gene2", "cor")
  cor_long <- cor_long %>% filter(!is.na(cor), gene1 < gene2)
  cor_long <- cor_long %>% mutate(abs_cor = abs(cor))
  cor_filtered <- cor_long %>% filter(abs_cor >= 0.7)
  if (nrow(cor_filtered) == 0) {
    cor_filtered <- cor_long %>% slice_max(abs_cor, n = 200, with_ties = FALSE)
  }
  write_csv(cor_filtered, file.path(res_dir, "network_edges_top50.csv"))
  degree_tbl <- cor_filtered %>%
    select(gene1, gene2) %>%
    pivot_longer(everything(), values_to = "gene") %>%
    count(gene, name = "degree") %>%
    arrange(desc(degree))
  write_csv(degree_tbl, file.path(res_dir, "network_hubs_top50.csv"))
  cor_plot_df <- cor_mat[genes_for_net, genes_for_net]
  cor_plot_df[lower.tri(cor_plot_df, diag = TRUE)] <- NA
  cor_plot_df <- cor_plot_df %>%
    as_tibble(rownames = "gene1") %>%
    pivot_longer(-gene1, names_to = "gene2", values_to = "cor") %>%
    filter(!is.na(cor))
  p_cor_heat <- ggplot(cor_plot_df, aes(gene1, gene2, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), oob = squish) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Top DE genes correlation (upper triangle)", x = "", y = "")
  ggsave(file.path(res_dir, "network_top50_heatmap.png"), p_cor_heat, width = 8, height = 7, dpi = 200)
} else {
  message("Not enough genes for network plot.")
}

# PCA after optional ComBat-seq (for visualization only)
if (!is.null(combat_counts)) {
  combat_log <- log2(combat_counts + 1)
  combat_top <- head(order(rowVars(combat_log, na.rm = TRUE), decreasing = TRUE), 2000)
  pca_bc <- prcomp(t(combat_log[combat_top, ]), center = TRUE, scale. = TRUE)
  pca_bc_df <- tibble(
    sample_id = colnames(combat_log),
    dataset = c(rep("GTEx", ncol(gtex)), rep("TCGA", ncol(tcga))),
    group = group,
    PC1 = pca_bc$x[, 1],
    PC2 = pca_bc$x[, 2]
  )
  p_bc <- ggplot(pca_bc_df, aes(PC1, PC2, color = dataset, shape = group)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(GTEx = "#1b9e77", TCGA = "#d95f02")) +
    guides(shape = guide_legend(title = "Sample type"))
  ggsave(file.path(res_dir, "pca_combatseq.png"), p_bc, width = 7, height = 5, dpi = 200)
}

# Survival analysis: stage + age
surv_complete <- survival_df %>%
  filter(!is.na(time_days), time_days > 0, !is.na(tumor_stage), tumor_stage != "") %>%
  mutate(tumor_stage = factor(tumor_stage))

cox_stage_age <- coxph(Surv(time_days, event) ~ tumor_stage + age_at_diagnosis, data = surv_complete)
cox_stage_sum <- summary(cox_stage_age)
cox_stage_tbl <- tibble(
  term = rownames(cox_stage_sum$coefficients),
  coef = cox_stage_sum$coefficients[, "coef"],
  HR = cox_stage_sum$conf.int[, "exp(coef)"],
  lower95 = cox_stage_sum$conf.int[, "lower .95"],
  upper95 = cox_stage_sum$conf.int[, "upper .95"],
  pvalue = cox_stage_sum$coefficients[, "Pr(>|z|)"]
)
write_csv(cox_stage_tbl, file.path(res_dir, "cox_stage_age.csv"))

km_stage_df <- km_step_df(surv_complete, "time_days", "event", "tumor_stage")
p_km_stage <- ggplot(km_stage_df, aes(time, surv, color = group)) +
  geom_step() +
  labs(x = "Days", y = "Survival probability", color = "Stage") +
  theme(legend.position = "right")
ggsave(file.path(res_dir, "km_by_stage.png"), p_km_stage, width = 7, height = 5, dpi = 200)

# Survival + single gene analysis for top N DE genes
top_genes <- head(de_tbl$gene_id, 5)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

cox_gene_all <- list()
single_stats_all <- list()

for (single_gene in top_genes) {
  if (!(single_gene %in% rownames(gtex_logtpm_common) && single_gene %in% rownames(tcga_logtpm_common))) next

  single_df <- bind_rows(
    tibble(
      dataset = "GTEx",
      sample_id = colnames(gtex_logtpm_common),
      group = "normal",
      expr = gtex_logtpm_common[single_gene, ]
    ),
    tibble(
      dataset = "TCGA",
      sample_id = colnames(tcga_logtpm_common),
      group = "tumor",
      expr = tcga_logtpm_common[single_gene, ],
      tumor_stage = colData(tcga)$tcga.gdc_cases.diagnoses.tumor_stage
    )
  )

  p_box <- ggplot(single_df, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = c(normal = "#1b9e77", tumor = "#d95f02")) +
    labs(title = paste0(single_gene, " expression"), x = "", y = "logTPM") +
    guides(fill = "none")
  ggsave(file.path(res_dir, paste0("single_gene_box_", safe_name(single_gene), ".png")), p_box, width = 5, height = 4, dpi = 200)

  tcga_stage_df <- single_df %>%
    filter(dataset == "TCGA", !is.na(tumor_stage), tumor_stage != "")
  if (nrow(tcga_stage_df) > 0) {
    p_stage <- ggplot(tcga_stage_df, aes(tumor_stage, expr, fill = tumor_stage)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
      labs(title = paste0(single_gene, " in TCGA by stage"), x = "Stage", y = "logTPM") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(res_dir, paste0("single_gene_stage_", safe_name(single_gene), ".png")), p_stage, width = 7, height = 4.5, dpi = 200)
  }

  stats_df <- single_df %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(expr, na.rm = TRUE),
      median = median(expr, na.rm = TRUE),
      sd = sd(expr, na.rm = TRUE),
      gene = single_gene,
      .groups = "drop"
    )
  single_stats_all[[single_gene]] <- stats_df

  surv_gene_df <- survival_df %>%
    left_join(single_df %>% select(sample_id, expr), by = "sample_id") %>%
    filter(!is.na(expr), !is.na(time_days), time_days > 0) %>%
    mutate(expr_group = if_else(expr >= median(expr), "High", "Low"))

  if (nrow(surv_gene_df) >= 10) {
    cox_gene <- coxph(Surv(time_days, event) ~ expr + age_at_diagnosis, data = surv_gene_df)
    cox_gene_sum <- summary(cox_gene)
    cox_gene_tbl <- tibble(
      gene = single_gene,
      term = rownames(cox_gene_sum$coefficients),
      coef = cox_gene_sum$coefficients[, "coef"],
      HR = cox_gene_sum$conf.int[, "exp(coef)"],
      lower95 = cox_gene_sum$conf.int[, "lower .95"],
      upper95 = cox_gene_sum$conf.int[, "upper .95"],
      pvalue = cox_gene_sum$coefficients[, "Pr(>|z|)"]
    )
    cox_gene_all[[single_gene]] <- cox_gene_tbl

    km_gene_df <- km_step_df(surv_gene_df, "time_days", "event", "expr_group")
    p_km_gene <- ggplot(km_gene_df, aes(time, surv, color = group)) +
      geom_step() +
      labs(x = "Days", y = "Survival probability", color = paste0(single_gene, " expression")) +
      theme(legend.position = "right")
    ggsave(file.path(res_dir, paste0("km_gene_", safe_name(single_gene), ".png")), p_km_gene, width = 7, height = 5, dpi = 200)
  }
}

if (length(cox_gene_all) > 0) {
  bind_rows(cox_gene_all) %>% write_csv(file.path(res_dir, "cox_top_genes.csv"))
}
if (length(single_stats_all) > 0) {
  bind_rows(single_stats_all) %>% write_csv(file.path(res_dir, "single_gene_stats.csv"))
}

message("EDA + DE + survival complete. Tables and plots are in the results/ directory.")
