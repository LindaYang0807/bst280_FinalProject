# BST280 Final Project: TCGA-BRCA ER+ vs ER- network analysis

This repository contains a full RNA-seq workflow to compare estrogen receptor-positive (ERpos) versus estrogen receptor-negative (ERneg) TCGA BRCA tumors and to characterize subtype-specific transcriptional networks.

## Data
- Input: `data/tcga_brca.rds` (expected as a `SummarizedExperiment`, or a list containing a count matrix plus clinical metadata). The metadata must include `tcga.xml_breast_carcinoma_estrogen_receptor_status` to classify samples as ERpos/ERneg.
- Output directories are created under `results/er_analysis/figures` and `results/er_analysis/tables`.

## How to run
From the repository root:
```bash
Rscript er_analysis_updated.R
```
- The script attempts to install required packages (`DESeq2`, `fgsea`, `msigdbr`, `org.Hs.eg.db`, `ggplot2`, `Rtsne`, etc.) and optional ones (`netZooR`, `apeglm`, `ashr`, `clusterProfiler`) if missing. First-time execution requires internet access for installation.
- PANDA uses `netZooR` when available; otherwise it falls back to correlation-based weights.
- To regenerate the narrative report, knit `Final_Project-copy.Rmd` (same pipeline with prose).

## Pipeline overview (er_analysis_updated.R)
1) Load counts/metadata, align shared samples, and filter for clear ER status (Positive/Negative), recoding to `ERneg`/`ERpos`.
2) Gene QC: strip Ensembl versions, drop duplicates, require counts >=10 in >=20% of samples, and remove genes lacking HGNC symbols.
3) DESeq2 analysis: run with parallel workers (no outlier replacement), obtain VST-normalized matrix with symbol mapping.
4) EDA: PCA on all genes and on top variable subsets (2k/500); t-SNE on top variable genes (2k/500).
5) Differential expression (ERpos vs ERneg): LFC shrinkage via `apeglm` or `ashr` when available; significance requires `padj < 0.05` and `|log2FC| >= 1`. Exports diagnostics (MA plot, volcano, p-value histogram) and summary metrics.
6) GSEA: `fgseaMultilevel` on GO Biological Process; saves full results, top 15 table, barplot, enrichment curves, leading-edge cilium genes, and per-sample cilium scores.
7) PANDA networks: build TF-target prior from MSigDB C3 TFT sets (limited to expressed genes, capped at 5k genes), infer ERpos and ERneg networks via `netZooR` or correlation fallback, and plot the top-weighted edges for each subtype.
8) Differential targeting: compare ERpos vs ERneg edge weights, export TF-level deltas, top TF target sets, barplot of top 20 TFs, and a heatmap of top TF-target deltas.

## Key outputs
- Figures (`results/er_analysis/figures`): PCA (`pca_vst_all.pdf`, `pca_vst_top2000.pdf`, `pca_vst_top500.pdf`), t-SNE (`tsne_tsne_top2000.pdf`, `tsne_tsne_top500.pdf`), DE plots (`volcano_erpos_vs_erneg.pdf`, `ma_plot_erpos_vs_erneg.pdf`, `pvalue_histogram.png`), GSEA (`fgsea_go_bp_top15.pdf`, `fgsea_go_bp_top15_plotEnrichment.pdf`), cilium score boxplot (`cilium_signature_score_boxplot.pdf`), PANDA networks (`panda_networks_top_edges.pdf`), differential targeting (`tf_differential_targeting_top20_bar.pdf`, `delta_topTFs_topTargets_heatmap.pdf`).
- Tables/objects (`results/er_analysis/tables`): symbol map (`ensembl_to_symbol_map.tsv`), PCA/t-SNE coordinates, DE results (`deseq2_erpos_vs_erneg_full_shrunk.csv`, `deseq2_summary_metrics.csv`, `deseq2_diagnostics.txt`, `ranked_genes_stat.tsv`, DESeq2 objects), GSEA outputs (`fgsea_go_bp_erpos_vs_erneg.csv`, `fgsea_go_bp_top15_labeled.csv`, leading-edge/cilium files, cilium scores), PANDA networks (`panda_network_erpos.rds`, `panda_network_erneg.rds`), differential targeting summaries (`differential_targeting_tf.csv`, `top_tf_target_sets.rds`).
