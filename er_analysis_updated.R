#!/usr/bin/env Rscript
# TCGA-BRCA ER+ vs ER− RNA-seq analysis
set.seed(20250219)
options(stringsAsFactors = FALSE)

# -------------------- Setup --------------------
required_pkgs <- c(
  "DESeq2", "ggplot2", "ggrepel", "dplyr", "tibble", "readr", "tidyr",
  "pheatmap", "fgsea", "msigdbr", "org.Hs.eg.db", "AnnotationDbi", "dplyr",
  "data.table", "BiocParallel", "matrixStats"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Package ", p, " not found; attempting install (requires internet).")
      tryCatch({
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(p, ask = FALSE, update = FALSE)
      }, error = function(e) {
        stop("Package ", p, " is required but could not be installed automatically. ",
             "Install it manually (internet needed) and re-run. Details: ", e$message)
      })
    }
    library(p, character.only = TRUE)
  }
}
install_if_missing(required_pkgs)

has_netZooR <- requireNamespace("netZooR", quietly = TRUE)
if (!has_netZooR) {
  message("netZooR not installed; PANDA will fall back to correlation-based weights.")
}
has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)
if (!has_clusterProfiler) {
  message("clusterProfiler not installed; TF-target enrichment will be skipped. Install manually to enable.")
}

workers <- max(1, parallel::detectCores() - 1)
bp_param <- BiocParallel::MulticoreParam(workers = workers)
use_parallel <- workers > 1
fgsea_nperm <- 2000  # reduce permutations for speed
max_network_genes <- 5000  # cap genes passed into PANDA/correlation
message("Parallel workers for DE/network: ", workers)

er_dir <- file.path("results", "er_analysis")
er_fig_dir <- file.path(er_dir, "figures")
er_tbl_dir <- file.path(er_dir, "tables")
dir.create(er_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(er_tbl_dir, recursive = TRUE, showWarnings = FALSE)
pca_log <- file.path(er_tbl_dir, "pca_diagnostics.txt")
if (file.exists(pca_log)) file.remove(pca_log)

# -------------------- Helpers --------------------
strip_version <- function(x) sub("\\..*$", "", x)
clean_er <- function(x) {
  y <- toupper(trimws(x))
  y <- ifelse(grepl("POS", y), "ERpos",
              ifelse(grepl("NEG", y), "ERneg", NA))
  y
}
write_placeholder_pdf <- function(path, message_line, width = 6, height = 4) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.text(message_line)
  grDevices::dev.off()
}

# -------------------- Load data --------------------
obj <- readRDS("tcga_brca.rds")

if (inherits(obj, "SummarizedExperiment")) {
  expr_mat <- as.matrix(SummarizedExperiment::assay(obj))
  meta <- as.data.frame(SummarizedExperiment::colData(obj))
} else if (is.list(obj)) {
  candidate_expr <- obj[intersect(names(obj), c("expression", "expr", "counts", "assay", "data"))]
  candidate_meta <- obj[intersect(names(obj), c("colData", "metadata", "clinical", "pheno", "samples"))]
  expr_mat <- if (length(candidate_expr) > 0) as.matrix(candidate_expr[[1]]) else stop("No expression matrix found.")
  meta <- if (length(candidate_meta) > 0) as.data.frame(candidate_meta[[1]]) else stop("No metadata found.")
} else if (is.matrix(obj) || is.data.frame(obj)) {
  expr_mat <- as.matrix(obj)
  meta <- stop("Metadata missing; cannot proceed without clinical data.")
} else {
  stop("Unsupported data structure in RDS.")
}

# Align samples
common_samples <- intersect(colnames(expr_mat), rownames(meta))
if (length(common_samples) < 4) stop("Too few overlapping samples.")
expr_mat <- expr_mat[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# -------------------- ER status filtering --------------------

er_col <- "tcga.xml_breast_carcinoma_estrogen_receptor_status"

if (!er_col %in% colnames(meta)) {
  stop(paste0("ER status column '", er_col, "' not found in metadata."))
}

# Inspect raw ER labels (useful for sanity check)
cat("Raw ER status distribution:\n")
print(table(meta[[er_col]], useNA = "ifany"))

# Keep only clear ER-positive / ER-negative samples; drop NA labels completely
keep <- !is.na(meta[[er_col]]) & meta[[er_col]] %in% c("Positive", "Negative")
if (sum(keep) < 2) {
  stop("Fewer than two samples have clear ER status (Positive/Negative).")
}

expr_mat <- expr_mat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

# Recode ER status
meta$ER_status <- ifelse(
  meta[[er_col]] == "Positive",
  "ERpos",
  "ERneg"
)

meta$ER_status <- factor(meta$ER_status, levels = c("ERneg", "ERpos"))
if (anyNA(meta$ER_status)) {
  stop("Detected NA ER_status entries after filtering; check metadata format.")
}
meta$ER_status <- droplevels(meta$ER_status)
if (!setequal(levels(meta$ER_status), c("ERneg", "ERpos"))) {
  stop("After filtering, ER_status must contain both ERneg and ERpos samples.")
}

cat("Samples retained after ER filter:", ncol(expr_mat), "\n")


# -------------------- Gene filtering --------------------
rownames(expr_mat) <- strip_version(rownames(expr_mat))
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]
keep_genes <- rowSums(expr_mat >= 10) >= ceiling(0.2 * ncol(expr_mat))
expr_filt <- expr_mat[keep_genes, ]
cat("Genes retained after low-expression filter:", nrow(expr_filt), "\n")

# Remove genes lacking a SYMBOL annotation
symbol_map_expr <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(expr_filt),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)
symbol_map_expr <- symbol_map_expr[!duplicated(symbol_map_expr$ENSEMBL), , drop = FALSE]
valid_symbol_ids <- symbol_map_expr$ENSEMBL[!is.na(symbol_map_expr$SYMBOL) & symbol_map_expr$SYMBOL != ""]
expr_filt <- expr_filt[rownames(expr_filt) %in% valid_symbol_ids, , drop = FALSE]
cat("Genes retained after removing entries without SYMBOL:", nrow(expr_filt), "\n")
if (nrow(expr_filt) == 0) {
  stop("No genes remain after filtering for SYMBOL annotations.")
}

# -------------------- DESeq2 setup & normalization --------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = expr_filt, colData = meta, design = ~ ER_status)

# More conservative fitting: don’t replace outliers (avoids artificially huge effects)
dds <- DESeq2::DESeq(
  dds,
  parallel = use_parallel,
  BPPARAM = bp_param,
  minReplicatesForReplace = Inf  # NEW
)

vsd <- DESeq2::vst(dds, blind = FALSE)
norm_mat <- SummarizedExperiment::assay(vsd)

# Map normalized matrix rownames (Ensembl) to gene symbols for plotting and clustering
map_norm <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(norm_mat),
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)
map_norm <- map_norm[!duplicated(map_norm$ENSEMBL), , drop = FALSE]
sym_vec <- setNames(map_norm$SYMBOL, map_norm$ENSEMBL)
sym_lookup <- sym_vec[rownames(norm_mat)]
num_symbol_na <- sum(is.na(sym_lookup) | sym_lookup == "")
cat("Rows lacking a mapped symbol:", num_symbol_na, "\n")
has_symbol <- !is.na(sym_lookup) & sym_lookup != ""
if (sum(has_symbol) == 0) {
  stop("No genes with valid SYMBOL IDs remain after filtering.")
}
if (sum(!has_symbol) > 0) {
  message("Filtering out ", sum(!has_symbol), " genes lacking SYMBOL annotations from normalized matrix.")
}
norm_mat_sym <- norm_mat[has_symbol, , drop = FALSE]
rownames(norm_mat_sym) <- sym_lookup[has_symbol]
rownames_norm_sym <- rownames(norm_mat_sym)
symbol_map_df <- data.frame(
  ensembl = rownames(norm_mat)[has_symbol],
  symbol = sym_lookup[has_symbol],
  stringsAsFactors = FALSE
)
readr::write_tsv(symbol_map_df, file.path(er_tbl_dir, "ensembl_to_symbol_map.tsv"))

# -------------------- PCA and clustering (all genes and variable-gene subsets) --------------------

# Helper: compute PCA, save coordinates, and create clearer ggplot output
make_pca_plot <- function(mat, meta, prefix = "pca", top_n = NA) {
  if (!is.na(top_n)) {
    if (nrow(mat) <= top_n) {
      mat_use <- mat
      sel_name <- paste0(prefix, "_top", nrow(mat))
    } else {
      v <- matrixStats::rowVars(mat, na.rm = TRUE)
      ord <- order(v, decreasing = TRUE)
      idx <- ord[seq_len(min(top_n, length(v)))]
      mat_use <- mat[idx, , drop = FALSE]
      sel_name <- paste0(prefix, "_top", top_n)
    }
  } else {
    mat_use <- mat
    sel_name <- paste0(prefix, "_all")
  }

  if (nrow(mat_use) < 3 || ncol(mat_use) < 3) {
    message("Skipping PCA for ", sel_name, ": too few genes/samples.")
    return(invisible(NULL))
  }

  pca_res <- tryCatch(prcomp(t(mat_use), scale. = TRUE), error = function(e) { message("PCA error: ", e$message); NULL })
  if (is.null(pca_res)) return(invisible(NULL))

  var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], ER_status = meta$ER_status, sample = rownames(pca_res$x), stringsAsFactors = FALSE)

  # Save coordinates for later inspection
  readr::write_tsv(pca_df, file.path(er_tbl_dir, paste0("pca_coords_", sel_name, ".tsv")))

  p <- ggplot(pca_df, aes(PC1, PC2, color = ER_status, label = sample)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(check_overlap = TRUE, vjust = 1.5, size = 2) +
    stat_ellipse(aes(color = ER_status), level = 0.68, linetype = 2, size = 0.6, show.legend = FALSE) +
    theme_bw() +
    labs(title = paste0("PCA (", sel_name, ")"),
         x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"))

  ggsave(file.path(er_fig_dir, paste0("pca_", sel_name, ".pdf")), p, width = 6, height = 5)
  invisible(p)
}

# Produce PCA plots for: all genes, top 2000 most variable, and top 500 most variable genes
make_pca_plot(norm_mat_sym, meta, prefix = "vst", top_n = NA)
make_pca_plot(norm_mat_sym, meta, prefix = "vst", top_n = 2000)
make_pca_plot(norm_mat_sym, meta, prefix = "vst", top_n = 500)


# -------------------- Differential expression (with shrinkage + effect-size filter) --------------------

res <- DESeq2::results(
  dds,
  contrast = c("ER_status", "ERpos", "ERneg"),
  independentFiltering = TRUE,
  alpha = 0.05
)

# Prefer shrunken LFCs: avoids tiny p but absurd LFCs for low-count genes
if (requireNamespace("apeglm", quietly = TRUE)) {
  message("Using apeglm for LFC shrinkage.")
  res_shrunk <- DESeq2::lfcShrink(
    dds,
    coef = "ER_status_ERpos_vs_ERneg",  # name from resultsNames(dds)
    type = "apeglm"
  )
  res_use <- res_shrunk
} else if (requireNamespace("ashr", quietly = TRUE)) {
  message("Using ashr for LFC shrinkage.")
  res_shrunk <- DESeq2::lfcShrink(
    dds,
    coef = "ER_status_ERpos_vs_ERneg",
    type = "ashr"
  )
  res_use <- res_shrunk
} else {
  message("No shrinkage package installed (apeglm / ashr); using raw DESeq2 results.")
  res_use <- res
}

res_tbl <- as.data.frame(res_use) %>%
  tibble::rownames_to_column("ensembl") %>%
  mutate(
    ensembl = strip_version(ensembl)
  ) %>%
  # drop rows with NA p-values
  filter(!is.na(pvalue)) %>%
  # conservative family-wise correction if you still want it
  mutate(padj_bonf = p.adjust(pvalue, method = "bonferroni"))

# Add a *biological* effect-size filter: e.g. |log2FC| ≥ 1 (2-fold change)
lfc_cutoff <- 1  # you can change to log2(1.5) if you want milder (≈1.5x) effects
res_tbl <- res_tbl %>%
  mutate(
    is_sig = (padj < 0.05) & (abs(log2FoldChange) >= lfc_cutoff)
  )

cat("Number of DE genes with padj < 0.05 and |log2FC| ≥ ", lfc_cutoff, ": ",
    sum(res_tbl$is_sig, na.rm = TRUE), "\n")

# Map Ensembl IDs to gene symbols; fall back to Ensembl ID when no symbol found
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_tbl$ensembl,
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

map_df <- map_df[!duplicated(map_df$ENSEMBL), , drop = FALSE]

res_tbl <- res_tbl %>%
  left_join(map_df, by = c("ensembl" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  mutate(gene = SYMBOL) %>%
  arrange(padj)  # sort by FDR

readr::write_csv(res_tbl, file.path(er_tbl_dir, "deseq2_erpos_vs_erneg_full_shrunk.csv"))

# Some shrinkers (e.g., apeglm) drop the Wald statistic; rebuild it for ranking
if (!"stat" %in% colnames(res_tbl)) {
  res_tbl <- res_tbl %>%
    mutate(
      stat = dplyr::if_else(
        !is.na(lfcSE) & lfcSE > 0,
        log2FoldChange / lfcSE,
        NA_real_
      )
    )
} else if (!is.numeric(res_tbl$stat)) {
  res_tbl$stat <- as.numeric(res_tbl$stat)
}


# Save DESeq2 objects for reproducibility
saveRDS(dds, file.path(er_tbl_dir, "dds_er_analysis.rds"))
saveRDS(res, file.path(er_tbl_dir, "deseq2_results_object.rds"))
saveRDS(vsd, file.path(er_tbl_dir, "vsd_er_analysis.rds"))
saveRDS(norm_mat, file.path(er_tbl_dir, "vst_normalized_matrix.rds"))

# -------------------- Basic sanity checks & diagnostic plots --------------------
n_genes_total <- nrow(res_tbl)
# Use Bonferroni-adjusted p-value for significance
n_sig <- sum(!is.na(res_tbl$padj_bonf) & res_tbl$padj_bonf < 0.05)
n_up <- sum(!is.na(res_tbl$padj_bonf) & res_tbl$padj_bonf < 0.05 & res_tbl$log2FoldChange > 1)
n_down <- sum(!is.na(res_tbl$padj_bonf) & res_tbl$padj_bonf < 0.05 & res_tbl$log2FoldChange < -1)
top_genes <- res_tbl %>% slice_head(n = 10) %>% pull(gene)
esr1_hit <- FALSE
if ("ESR1" %in% res_tbl$gene) {
  esr1_row <- res_tbl %>% filter(gene == "ESR1")
  esr1_hit <- ifelse(nrow(esr1_row) > 0 && !is.na(esr1_row$padj_bonf) && esr1_row$padj_bonf < 0.05, TRUE, FALSE)
}

summary_tbl <- tibble::tibble(
  metric = c("n_genes_total", "n_sig", "n_up", "n_down", "esr1_significant"),
  value = c(n_genes_total, n_sig, n_up, n_down, as.character(esr1_hit))
)
readr::write_csv(summary_tbl, file.path(er_tbl_dir, "deseq2_summary_metrics.csv"))

# MA plot (saved as PDF)
try({
  pdf(file.path(er_fig_dir, "ma_plot_erpos_vs_erneg.pdf"), width = 6, height = 5)
  DESeq2::plotMA(res, main = "MA plot: ERpos vs ERneg", ylim = c(-5, 5))
  dev.off()
}, silent = TRUE)

# P-value distribution histogram (saved as PNG)
try({
  pval_vec <- res$pvalue
  pval_vec <- pval_vec[!is.na(pval_vec)]
  if (length(pval_vec) > 0) {
    png(file.path(er_fig_dir, "pvalue_histogram.png"), width = 800, height = 600)
    hist(pval_vec, breaks = 50, main = "P-value distribution", xlab = "p-value", col = "grey80")
    dev.off()
  }
}, silent = TRUE)

# Quick textual diagnostics
diag_lines <- c(
  paste0("DESeq2: total genes tested = ", n_genes_total),
  paste0("DESeq2: significant (padj_bonf < 0.05) = ", n_sig),
  paste0("DESeq2: upregulated in ERpos (padj_bonf < 0.05 & log2FC>1) = ", n_up),
  paste0("DESeq2: downregulated in ERpos (padj_bonf < 0.05 & log2FC<-1) = ", n_down),
  paste0("ESR1 significant (padj_bonf < 0.05) = ", esr1_hit),
  paste0("Top genes (by padj_bonf): ", paste(head(top_genes, 10), collapse = ", "))
)
readr::write_lines(diag_lines, file.path(er_tbl_dir, "deseq2_diagnostics.txt"))

# Create a deduplicated ranked gene list by gene symbol (choose the row with largest |stat| when symbols duplicate)
ranked_genes <- res_tbl %>%
  filter(!is.na(stat)) %>%
  group_by(gene) %>%
  slice_max(order_by = abs(stat), n = 1) %>%
  ungroup() %>%
  mutate(rank = rank(-stat, ties.method = "average")) %>%
  arrange(rank)
readr::write_tsv(ranked_genes[, c("gene", "stat")], file.path(er_tbl_dir, "ranked_genes_stat.tsv"))

volcano_p_thresh <- 0.05

volcano_data <- res_tbl %>%
  mutate(
    direction = case_when(
      !is.na(padj_bonf) & padj_bonf < volcano_p_thresh & log2FoldChange > 1 ~ "Up",
      !is.na(padj_bonf) & padj_bonf < volcano_p_thresh & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    neglog10_padj_bonf = ifelse(
      is.na(padj_bonf),
      NA_real_,
      -log10(pmax(padj_bonf, .Machine$double.xmin))
    )
  )
label_data <- volcano_data %>%
  filter(direction %in% c("Up", "Down")) %>%
  arrange(padj_bonf) %>%
  slice_head(n = 20)

p_volcano <- ggplot(volcano_data, aes(log2FoldChange, neglog10_padj_bonf)) +
  geom_point(aes(color = direction), alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = -log10(volcano_p_thresh), linetype = 2, color = "grey50", size = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, color = "grey70", size = 0.3) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 2.4,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(Up = "firebrick", Down = "royalblue3", NS = "grey80")) +
  labs(
    title = "ERpos vs ERneg DESeq2 (Bonferroni)",
    x = "log2 Fold Change (ERpos/ERneg)",
    y = "-log10(padj_bonf)",
    color = "Direction"
  ) +
  theme_bw()
ggsave(file.path(er_fig_dir, "volcano_erpos_vs_erneg.pdf"), p_volcano, width = 6, height = 5)


# -------------------- GSEA (GO BP) --------------------
gene_ranks <- ranked_genes %>%
  dplyr::select(gene, stat) %>%
  tibble::deframe()
gene_ranks <- gene_ranks[!is.na(gene_ranks)]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

go_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
) %>%
  dplyr::select(gs_name, gene_symbol)
go_list <- split(go_sets$gene_symbol, go_sets$gs_name)
fgsea_res <- fgsea::fgsea(pathways = go_list, stats = gene_ranks, nperm = fgsea_nperm)
fgsea_tbl <- fgsea_res %>%
  arrange(padj) %>%
  as_tibble()
readr::write_csv(fgsea_tbl, file.path(er_tbl_dir, "fgsea_go_bp_erpos_vs_erneg.csv"))
fgsea_top15 <- fgsea_tbl %>%
  filter(!is.na(padj)) %>%
  slice_head(n = 15) %>%
  mutate(
    rank = dplyr::row_number(),
    pathway_label = paste0(rank, ". ", pathway)
  )
if (nrow(fgsea_top15) > 0) {
  readr::write_csv(fgsea_top15, file.path(er_tbl_dir, "fgsea_go_bp_top15_labeled.csv"))
  fgsea_top15_plot <- fgsea_top15 %>%
    mutate(pathway_label = factor(pathway_label, levels = rev(pathway_label)))
  p_fgsea_top15 <- ggplot(fgsea_top15_plot, aes(x = pathway_label, y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "royalblue3", mid = "grey90", high = "firebrick", midpoint = 0) +
    labs(
      title = "Top 15 GO BP pathways (fgsea)",
      x = NULL,
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 8)
    )
  ggsave(
    filename = file.path(er_fig_dir, "fgsea_go_bp_top15.pdf"),
    plot = p_fgsea_top15,
    width = 7,
    height = 5,
    useDingbats = FALSE
  )
  enrich_pdf <- file.path(er_fig_dir, "fgsea_go_bp_top15_plotEnrichment.pdf")
  pdf_open <- FALSE
  tryCatch({
    grDevices::pdf(enrich_pdf, width = 6, height = 4, useDingbats = FALSE)
    pdf_open <- TRUE
    for (i in seq_len(nrow(fgsea_top15))) {
      path <- fgsea_top15$pathway[i]
      if (!path %in% names(go_list)) next
      label <- fgsea_top15$pathway_label[i]
      p_enrich <- fgsea::plotEnrichment(go_list[[path]], gene_ranks) +
        ggtitle(label)
      print(p_enrich)
    }
    grDevices::dev.off()
    pdf_open <- FALSE
  }, error = function(e) {
    if (pdf_open && grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    warning("Failed to save fgsea plotEnrichment PDF: ", e$message)
  })
}


# -------------------- PANDA network inference --------------------
# Build motif prior from TF target gene sets (MSigDB C3:TFT) limited to expressed genes.
tft_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "C3",
  subcategory = "TFT:GTRD"
) %>%
  mutate(tf = gsub("_.*", "", gs_name)) %>%
  filter(gene_symbol %in% rownames_norm_sym)

motif_prior <- tft_sets %>%
  dplyr::select(tf, target = gene_symbol) %>%
  distinct() %>%
  mutate(score = 1)

ppi_prior <- motif_prior %>%
  distinct(tf) %>%
  transmute(protein1 = tf, protein2 = tf, score = 1)

expr_by_status <- list(
  ERpos = norm_mat_sym[, meta$ER_status == "ERpos", drop = FALSE],
  ERneg = norm_mat_sym[, meta$ER_status == "ERneg", drop = FALSE]
)

panda_step_records <- list()
panda_networks <- list()

for (status_label in names(expr_by_status)) {
  message("----- PANDA steps for ", status_label, " -----")
  expr_sub <- expr_by_status[[status_label]]

  # Step 1: Harmonize gene universe
  shared_genes <- intersect(rownames(expr_sub), motif_prior$target)
  if (length(shared_genes) > max_network_genes) {
    var_order <- order(
      matrixStats::rowVars(expr_sub[shared_genes, , drop = FALSE], na.rm = TRUE),
      decreasing = TRUE
    )
    shared_genes <- shared_genes[var_order[seq_len(max_network_genes)]]
  }
  message("Step 1: ", length(shared_genes), " shared genes retained.")

  # Step 2: Subset expression and motif prior
  expr_use <- expr_sub[shared_genes, , drop = FALSE]
  motif_use <- motif_prior %>% filter(target %in% shared_genes)
  message("Step 2: motif subset has ", nrow(motif_use), " TF-target pairs.")

  net_df <- tibble::tibble(TF = character(), Gene = character(), Weight = numeric())
  method_used <- "not_run"
  tf_in_expr <- NULL
  targets_in_expr <- NULL

  if (length(shared_genes) == 0 || nrow(motif_use) == 0) {
    message("No overlapping genes or motif entries for ", status_label, "; skipping PANDA.")
    method_used <- "no_motif_overlap"
  } else if (has_netZooR) {
    # Step 3A: Run netZooR PANDA
    message("Step 3A: running netZooR PANDA for ", status_label)
    motif_netzoo <- motif_use
    colnames(motif_netzoo) <- c("TF", "Gene", "Strength")
    ppi_netzoo <- ppi_prior
    colnames(ppi_netzoo) <- c("Protein1", "Protein2", "Weight")
    panda_net <- netZooR::panda(
      motif = motif_netzoo,
      ppi = ppi_netzoo,
      expression = expr_use
    )
    net_df <- netZooR::pandaToDataFrame(panda_net)
    colnames(net_df) <- c("TF", "Gene", "Weight")
    method_used <- "netZooR"
  } else {
    # Step 3B: Correlation / proxy fallback
    tf_in_expr      <- intersect(unique(motif_use$tf), rownames(expr_sub))
    targets_in_expr <- intersect(motif_use$target, rownames(expr_use))
    message("Step 3B: fallback with ", length(tf_in_expr), " TFs and ", length(targets_in_expr), " targets.")

    if (length(targets_in_expr) < 1) {
      message("Skipping PANDA fallback for ", status_label, ": no expressed targets.")
      method_used <- "insufficient_targets"
    } else if (length(tf_in_expr) >= 1) {
      cor_mat <- stats::cor(
        t(expr_use[targets_in_expr, , drop = FALSE]),
        t(expr_sub[tf_in_expr,      , drop = FALSE]),
        use = "pairwise.complete.obs"
      )
      cor_df <- as.data.frame(as.table(cor_mat))
      if (ncol(cor_df) == 3) {
        colnames(cor_df) <- c("target", "tf", "weight")
        net_df <- cor_df %>%
          inner_join(motif_use, by = c("tf", "target")) %>%
          transmute(TF = tf, Gene = target, Weight = weight)
        method_used <- "correlation"
      } else {
        message("Unexpected correlation table shape; returning empty network for ", status_label)
        method_used <- "correlation_error"
      }
    } else {
      message("No TF genes detected for ", status_label, "; approximating from target sets.")
      tf_targets <- split(motif_use$target, motif_use$tf)
      tf_targets <- lapply(tf_targets, function(tg) intersect(tg, targets_in_expr))
      tf_targets <- tf_targets[lengths(tf_targets) > 0]
      if (length(tf_targets) == 0) {
        warning("No TF target sets passed filtering for ", status_label, ".")
        method_used <- "proxy_empty"
      } else {
        proxy_list <- lapply(names(tf_targets), function(tf) {
          tg <- tf_targets[[tf]]
          expr_block <- expr_use[tg, , drop = FALSE]
          if (nrow(expr_block) == 1) {
            weights <- 0
          } else {
            pseudo_tf <- colMeans(expr_block, na.rm = TRUE)
            weights <- apply(expr_block, 1, function(gexpr) {
              stats::cor(gexpr, pseudo_tf, use = "pairwise.complete.obs")
            })
          }
          tibble::tibble(TF = tf, Gene = tg, Weight = as.numeric(weights))
        })
        net_df <- dplyr::bind_rows(proxy_list) %>%
          dplyr::filter(!is.na(Weight))
        method_used <- "target_proxy"
      }
    }
  }

  panda_step_records[[status_label]] <- list(
    shared_genes = shared_genes,
    expr_subset = expr_use,
    motif_subset = motif_use,
    tf_in_expr = tf_in_expr,
    targets_in_expr = targets_in_expr,
    method = method_used,
    network = net_df
  )
  panda_networks[[status_label]] <- net_df
  message("Completed PANDA for ", status_label, " using method ", method_used,
          "; edges: ", nrow(net_df))
}

net_erpos <- panda_networks[["ERpos"]]
if (is.null(net_erpos)) {
  net_erpos <- tibble::tibble(TF = character(), Gene = character(), Weight = numeric())
}
net_erneg <- panda_networks[["ERneg"]]
if (is.null(net_erneg)) {
  net_erneg <- tibble::tibble(TF = character(), Gene = character(), Weight = numeric())
}

saveRDS(net_erpos, file.path(er_tbl_dir, "panda_network_erpos.rds"))
saveRDS(net_erneg, file.path(er_tbl_dir, "panda_network_erneg.rds"))


# -------------------- Differential targeting --------------------
net_merged <- net_erpos %>%
  rename(Weight_pos = Weight) %>%
  inner_join(net_erneg %>% rename(Weight_neg = Weight), 
             by = c("TF", "Gene"))%>%
  mutate(delta = Weight_pos - Weight_neg)

tf_diff <- net_merged %>%
  group_by(TF) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    n_targets = n()
  ) %>%
  arrange(desc(abs(mean_delta)))
readr::write_csv(tf_diff, file.path(er_tbl_dir, "differential_targeting_tf.csv")) 


tf_diff <- tf_diff %>%
  mutate(direction = ifelse(mean_delta > 0, "ERpos_higher", "ERneg_higher"))


# barplot of top 20 TFs by |mean_delta|
# Figure 1: Top TFs differential targeting barplot 
tf_top20 <- tf_diff %>%
  slice_max(order_by = abs(mean_delta), n = 20)

p_tf_bar <- ggplot(tf_top20,
                   aes(x = reorder(TF, mean_delta),
                       y = mean_delta,
                       fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top differentially targeted TFs (ERpos vs ERneg)",
    x     = "TF",
    y     = "Mean delta (Weight_pos - Weight_neg)",
    fill  = "Direction"
  )

ggsave(
  filename = file.path(er_fig_dir, "tf_differential_targeting_top20_bar.pdf"),
  plot     = p_tf_bar,
  width    = 7,
  height   = 5
)


# save top TFs
top_tfs <- tf_diff %>% slice_head(n = 15) %>% pull(TF)
tf_targets_list <- lapply(top_tfs, function(tf) {
  net_merged %>% filter(TF == tf) %>% arrange(desc(delta)) %>% pull(Gene)
})
names(tf_targets_list) <- top_tfs
saveRDS(tf_targets_list, file.path(er_tbl_dir, "top_tf_target_sets.rds"))


# Heatmap: For each top TF, keep top N targets by |delta|
# Figure 2: Heatmap of delta for top TFs and top targets 
top_n_targets_per_tf <- 50

delta_long <- net_merged %>%
  filter(TF %in% top_tfs) %>%
  group_by(TF) %>%
  slice_max(order_by = abs(delta), n = top_n_targets_per_tf, with_ties = FALSE) %>%
  ungroup()

# wide matrix: rows = TF, cols = Gene, values = delta
delta_wide <- delta_long %>%
  dplyr::select(TF, Gene, delta) %>%  # AnnotationDbi also defines select(); qualify to keep tibble method
  tidyr::pivot_wider(names_from = Gene, values_from = delta)

delta_mat <- delta_wide %>%
  as.data.frame()

rownames(delta_mat) <- delta_mat$TF
delta_mat$TF <- NULL
delta_mat <- as.matrix(delta_mat)
# Heatmap clustering cannot handle NAs; genes missing for a TF get weight 0
delta_mat[is.na(delta_mat)] <- 0

# Heatmap (pheatmap)
if (nrow(delta_mat) > 1 && ncol(delta_mat) > 1) {
  pheatmap::pheatmap(
    delta_mat,
    scale          = "none",
    clustering_method = "ward.D2",
    main           = "Differential targeting (delta) for top TFs and targets",
    filename       = file.path(er_fig_dir, "delta_topTFs_topTargets_heatmap.pdf"),
    width          = 10,
    height         = 7,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    fontsize_col   = 6
  )
}

