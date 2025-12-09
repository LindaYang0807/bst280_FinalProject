# BST280 Final Project

本项目提供 GTEx 乳腺（正常）与 TCGA BRCA（肿瘤）的探索性分析、差异表达与生存分析流水线。

## 内容概览
- 数据：`data/gtex_breast.rds`（GTEx 正常乳腺），`data/tcga_brca.rds`（TCGA BRCA 肿瘤），均为 `RangedSummarizedExperiment`，包含 `raw_counts` 与 `logtpm`。
- 脚本：`scripts/eda_gtex_tcga.R`，一键完成 EDA、差异表达、基础生存分析，结果输出到 `results/`。
- 依赖：SummarizedExperiment、tidyverse、matrixStats、patchwork、survival。若安装了 limma/edgeR/sva 将自动使用 limma-voom 与可选 ComBat_seq；否则退化为 log-CPM + Welch t 检验（已在当前运行中触发）。

运行：在仓库根目录执行 `Rscript scripts/eda_gtex_tcga.R`。

## 主要步骤（脚本内部）
1) 元数据汇总：性别、年龄、样本类型、分期、存活状态等，写入对应 CSV。
2) 基因 ID 去版本号并合并重复，保留 GTEx/TCGA 交集。
3) 质量分布：文库大小直方图、每样本中位数 logTPM 分布。
4) PCA：共同基因的前 2,000 高变基因，查看 GTEx 正常 vs TCGA 肿瘤分离（`results/pzned.png`）。
5) 差异表达：GTEx 正常 vs TCGA 肿瘤。
   - 若有 limma/edgeR：limma-voom；否则使用 log-CPM + Welch t（当前运行）。
   - 输出表 `results/de_limma_voom.csv`，火山图 `results/de_volcano.png`，MA 图 `results/de_MA.png`。
6) 生存分析：
   - 分期+年龄 Cox：`results/cox_stage_age.csv`；KM 分期曲线 `results/km_by_stage.png`。
   - 顶 DE 基因（当前为 ENSG00000000003）的高低表达分组 Cox：`results/cox_top_gene.csv`；KM 曲线 `results/km_top_gene.png`。
   - 生存输入表：`results/tcga_survival_ready.csv`。

## 已生成的结果（当前运行）
- 差异表达：`de_limma_voom.csv`（因缺少 limma/edgeR，使用 log-CPM + Welch t），`de_volcano.png`，`de_MA.png`。
- 分期生存：`cox_stage_age.csv`，`km_by_stage.png`。
- 顶基因生存：`cox_top_gene.csv`，`km_top_gene.png`。
- 元数据计数/EDA：`gtex_sex_counts.csv`，`gtex_age_counts.csv`，`tcga_sample_type_counts.csv`，`tcga_stage_counts.csv`，`tcga_vital_counts.csv`，`tcga_survival_ready.csv`，以及基础分布图（文库大小、logTPM 中位数、PCA）。

## 结果意义/后续可用性
- 差异表达结果提供肿瘤 vs 正常的候选基因列表，可用于功能注释、通路富集或后续生物标志物筛选。
- 生存分析给出分期与年龄的风险比，并用顶 DE 基因的表达高低分组展示潜在预后关联，可扩展到自定义基因或通路评分。
- 脚本是可重复的一键流程：补装 limma/edgeR/sva 即可切换到更标准的 voom/ComBat_seq 流程；可根据需求替换或增加基因签名、分组变量、统计模型。
