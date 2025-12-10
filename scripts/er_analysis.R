#!/usr/bin/env Rscript
# TCGA-BRCA ER+ vs ER− RNA-seq analysis
set.seed(20250219)
options(stringsAsFactors = FALSE)

# -------------------- Setup --------------------
required_pkgs <- c(
  "DESeq2", "ggplot2", "dplyr", "tibble", "readr", "tidyr",
  "pheatmap", "fgsea", "msigdbr", "org.Hs.eg.db", "AnnotationDbi",
  "data.table", "netZooR"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
    library(p, character.only = TRUE)
  }
}
install_if_missing(required_pkgs)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -------------------- Helpers --------------------
strip_version <- function(x) sub("\\..*$", "", x)
clean_er <- function(x) {
  y <- toupper(trimws(x))
  y <- ifelse(grepl("POS", y), "ERpos",
              ifelse(grepl("NEG", y), "ERneg", NA))
  y
}

# -------------------- Load data --------------------
obj <- readRDS("data/tcga_brca.rds")

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

# Keep only clear ER-positive / ER-negative samples
keep <- meta[[er_col]] %in% c("Positive", "Negative")

expr_mat <- expr_mat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

# Recode ER status
meta$ER_status <- ifelse(
  meta[[er_col]] == "Positive",
  "ERpos",
  "ERneg"
)

meta$ER_status <- factor(meta$ER_status, levels = c("ERneg", "ERpos"))

cat("Samples retained after ER filter:", ncol(expr_mat), "\n")


# -------------------- Gene filtering --------------------
rownames(expr_mat) <- strip_version(rownames(expr_mat))
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]
keep_genes <- rowSums(expr_mat >= 10) >= ceiling(0.2 * ncol(expr_mat))
expr_filt <- expr_mat[keep_genes, ]
cat("Genes retained after low-expression filter:", nrow(expr_filt), "\n")

# -------------------- DESeq2 setup & normalization --------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = expr_filt, colData = meta, design = ~ ER_status)
dds <- DESeq2::DESeq(dds)  # DESeq2 chosen for count-based RNA-seq with stable shrinkage
vsd <- DESeq2::vst(dds)
norm_mat <- SummarizedExperiment::assay(vsd)

# -------------------- PCA and clustering (all genes) --------------------
pca_res <- prcomp(t(norm_mat), scale. = TRUE)
pca_df <- data.frame(pca_res$x[, 1:2], ER_status = meta$ER_status, sample = rownames(pca_res$x))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = ER_status, label = sample)) +
  geom_point(size = 3) + geom_text(vjust = 1.5, size = 2) +
  theme_bw() + labs(title = "PCA on vst counts", color = "ER status")
ggsave("results/figures/pca_all_genes_er_status.pdf", p_pca, width = 6, height = 5)

km_res <- kmeans(t(norm_mat), centers = 2, nstart = 50)
km_table <- as.data.frame.matrix(table(km_res$cluster, meta$ER_status))
km_table$cluster <- rownames(km_table)
readr::write_csv(km_table, "results/tables/kmeans_cluster_vs_er.csv")

hc_res <- hclust(dist(t(norm_mat)), method = "ward.D2")
hc_clusters <- cutree(hc_res, k = 2)
hc_table <- as.data.frame.matrix(table(hc_clusters, meta$ER_status))
hc_table$cluster <- rownames(hc_table)
readr::write_csv(hc_table, "results/tables/hclust_cluster_vs_er.csv")
pdf("results/figures/hclust_dendrogram_all_genes.pdf", width = 6, height = 5)
plot(hc_res, labels = meta$ER_status, main = "Hierarchical clustering (all genes)")
abline(h = tail(hc_res$height, 1) / 2, col = "red", lty = 2)
dev.off()

# -------------------- Differential expression --------------------
res <- DESeq2::results(dds, contrast = c("ER_status", "ERpos", "ERneg"))
res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
readr::write_csv(res_tbl, "results/tables/deseq2_erpos_vs_erneg_full.csv")

ranked_genes <- res_tbl %>%
  filter(!is.na(stat)) %>%
  mutate(rank = rank(-stat, ties.method = "average")) %>%
  arrange(rank)
readr::write_tsv(ranked_genes[, c("gene", "stat")], "results/tables/ranked_genes_stat.tsv")

p_volcano <- res_tbl %>%
  mutate(sig = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, "Significant", "NS")) %>%
  ggplot(aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c(NS = "grey70", Significant = "firebrick")) +
  theme_bw() + labs(title = "ERpos vs ERneg DESeq2", color = "")
ggsave("results/figures/volcano_erpos_vs_erneg.pdf", p_volcano, width = 6, height = 5)

# -------------------- PCA & clustering on DE genes --------------------
sig_genes <- res_tbl %>%
  filter(!is.na(padj), padj < 0.05) %>%
  pull(gene)
if (length(sig_genes) >= 3) {
  norm_sig <- norm_mat[sig_genes, , drop = FALSE]
  pca_sig <- prcomp(t(norm_sig), scale. = TRUE)
  pca_sig_df <- data.frame(pca_sig$x[, 1:2], ER_status = meta$ER_status, sample = rownames(pca_sig$x))
  p_pca_sig <- ggplot(pca_sig_df, aes(PC1, PC2, color = ER_status, label = sample)) +
    geom_point(size = 3) + geom_text(vjust = 1.5, size = 2) +
    theme_bw() + labs(title = "PCA on significant DE genes", color = "ER status")
  ggsave("results/figures/pca_sig_genes.pdf", p_pca_sig, width = 6, height = 5)
  
  ann_df <- data.frame(ER_status = meta$ER_status)
  rownames(ann_df) <- rownames(meta)
  pheatmap(norm_sig, scale = "row", annotation_col = ann_df,
           show_rownames = FALSE, fontsize_col = 8,
           main = "DE genes heatmap (vst counts)",
           filename = "results/figures/heatmap_sig_genes.pdf")
  # Comment: Subtype separation should tighten here if DE genes capture ER biology.
}

# -------------------- GSEA (GO BP) --------------------
gene_symbols <- res_tbl$gene
gene_ranks <- res_tbl$stat
names(gene_ranks) <- gene_symbols
gene_ranks <- gene_ranks[!is.na(gene_ranks)]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

go_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  dplyr::select(gs_name, gene_symbol)
go_list <- split(go_sets$gene_symbol, go_sets$gs_name)
fgsea_res <- fgsea::fgsea(pathways = go_list, stats = gene_ranks, nperm = 5000)
fgsea_tbl <- fgsea_res %>%
  arrange(padj) %>%
  as_tibble()
readr::write_csv(fgsea_tbl, "results/tables/fgsea_go_bp_erpos_vs_erneg.csv")
top_path <- fgsea_tbl$pathway[1]
if (!is.null(top_path) && !is.na(top_path)) {
  pdf("results/figures/fgsea_enrichment_top_pathway.pdf", width = 6, height = 5)
  plotEnrichment(go_list[[top_path]], gene_ranks) + ggtitle(top_path)
  dev.off()
}

# -------------------- PANDA network inference --------------------
# Build motif prior from TF target gene sets (MSigDB C3:TFT) limited to expressed genes.
tft_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "C3",
  subcollection = "TFT:GTRD"
) %>%
  mutate(tf = gsub("_.*", "", gs_name)) %>%
  filter(gene_symbol %in% rownames(expr_filt))
motif_prior <- tft_sets %>%
  dplyr::select(tf, target = gene_symbol) %>%
  distinct() %>%
  mutate(score = 1)

ppi_prior <- motif_prior %>%
  distinct(tf) %>%
  transmute(protein1 = tf, protein2 = tf, score = 1)

panda_input <- function(expr_sub, motif_df, ppi_df) {
  shared_genes <- intersect(rownames(expr_sub), motif_df$target)
  expr_use <- expr_sub[shared_genes, , drop = FALSE]
  motif_use <- motif_df %>% filter(target %in% shared_genes)
  if (requireNamespace("netZooR", quietly = TRUE)) {
    colnames(motif_use) <- c("TF", "Gene", "Strength")
    colnames(ppi_df) <- c("Protein1", "Protein2", "Weight")
    panda_net <- netZooR::panda(motif = motif_use, ppi = ppi_df, expression = expr_use)
    net_df <- netZooR::pandaToDataFrame(panda_net)
  } else {
    # Fallback: PANDA-like weight = Pearson correlation within motif prior
    cor_mat <- cor(t(expr_use[motif_use$target, , drop = FALSE]),
                   t(expr_use[intersect(motif_use$tf, rownames(expr_use)), , drop = FALSE]),
                   use = "pairwise.complete.obs")
    cor_df <- as.data.frame(as.table(cor_mat))
    colnames(cor_df) <- c("target", "tf", "weight")
    net_df <- cor_df %>%
      inner_join(motif_use, by = c("tf", "target")) %>%
      transmute(TF = tf, Gene = target, Weight = weight)
  }
  net_df
}

expr_erpos <- norm_mat[, meta$ER_status == "ERpos", drop = FALSE]
expr_erneg <- norm_mat[, meta$ER_status == "ERneg", drop = FALSE]

net_erpos <- panda_input(expr_erpos, motif_prior, ppi_prior)
net_erneg <- panda_input(expr_erneg, motif_prior, ppi_prior)

saveRDS(net_erpos, "results/tables/panda_network_erpos.rds")
saveRDS(net_erneg, "results/tables/panda_network_erneg.rds")

# -------------------- Differential targeting --------------------
net_merged <- net_erpos %>%
  rename(Weight_pos = Weight) %>%
  inner_join(net_erneg %>% rename(Weight_neg = Weight), by = c("TF", "Gene"))
net_merged <- net_merged %>%
  mutate(delta = Weight_pos - Weight_neg)

tf_diff <- net_merged %>%
  group_by(TF) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    n_targets = n()
  ) %>%
  arrange(desc(abs(mean_delta)))
readr::write_csv(tf_diff, "results/tables/differential_targeting_tf.csv")

top_tfs <- tf_diff %>% slice_head(n = 5) %>% pull(TF)
tf_targets_list <- lapply(top_tfs, function(tf) {
  net_merged %>% filter(TF == tf) %>% arrange(desc(delta)) %>% pull(Gene)
})
names(tf_targets_list) <- top_tfs
saveRDS(tf_targets_list, "results/tables/top_tf_target_sets.rds")

# -------------------- GSEA/ORA on TF target sets --------------------
go_term2gene <- go_sets %>% dplyr::select(gs_name, gene_symbol)
tf_ora_results <- lapply(top_tfs, function(tf) {
  tg <- tf_targets_list[[tf]]
  enr <- clusterProfiler::enricher(
    gene = tg,
    TERM2GENE = go_term2gene,
    minGSSize = 10,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  if (!is.null(enr)) {
    res <- as.data.frame(enr@result)
    res$TF <- tf
    res
  } else {
    NULL
  }
})
tf_ora_tbl <- dplyr::bind_rows(tf_ora_results)
if (nrow(tf_ora_tbl) > 0) {
  readr::write_csv(tf_ora_tbl, "results/tables/tf_target_go_enrichment.csv")
}

# -------------------- Overlap comparison --------------------
if (nrow(tf_ora_tbl) > 0) {
  de_gsea_paths <- fgsea_tbl$pathway[fgsea_tbl$padj < 0.05]
  overlap <- tf_ora_tbl %>%
    filter(qvalue < 0.05) %>%
    mutate(in_de_gsea = gs_name %in% de_gsea_paths)
  readr::write_csv(overlap, "results/tables/tf_target_go_overlap_with_de_gsea.csv")
}

# -------------------- Gene deep dive (example: ESR1) --------------------
# ESR1 (Estrogen receptor 1) encodes a nuclear hormone receptor that mediates estrogen signaling;
# it is the canonical driver of ER+ breast cancer and often amplified or overexpressed in this subtype.
# ER+ tumors rely on ESR1 activity; mutations (e.g., Y537S, D538G) confer endocrine resistance.
# To look up ESR1 on NCBI: visit https://www.ncbi.nlm.nih.gov/gene/2099
# For UCSC Genome Browser track: open https://genome.ucsc.edu/, select hg38, search "ESR1",
# then add a custom track with BED coordinates of ESR1 enhancers or ATAC-seq peaks to visualize context.

sessionInfo()
