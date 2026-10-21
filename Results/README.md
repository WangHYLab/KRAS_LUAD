# Processed results and sample metadata

This directory contains processed tables used by the analyses described in the [main README](../README.md). Raw-data access and dataset-version information are documented there.

**Terminology:** CIIR refers to **cancer-cell-intrinsic immune-related genes in KRAS-mutant LUAD**. The integrated gene tables below identify CIIR genes and their immune associations; KRAS responsiveness is assessed separately using the KRAS-inhibitor perturbation results. PDE4D is described in the manuscript as a **KRAS-responsive CIIR candidate**.

The saved `negative`/`positive` gene labels denote negative/positive associations with immune-cell abundance. Some plotting sections use the corresponding labels `suppressive`/`promoting`; the table descriptions here use association terminology. Immune-hot/immune-cold subtypes and PDE4D-high/PDE4D-low groups refer to their respective grouping methods.

## Sample lists and grouping tables

| File | Source | Generating script / method |
| --- | --- | --- |
| [`sk.csv`](sk.csv) | TCGA-LUAD expression and mutation data | [`1.get_data.R`](../Scripts/1.get_data.R): match expression and mutation barcodes using their first 16 characters, apply the recorded duplicate-barcode selection, retain primary-tumour barcodes ending in `01A`, and assign KRAS status. Supplied table: 502 samples, with `KRAS=1` for 138 mutant and `KRAS=0` for 364 wild-type samples. |
| [`mutedata.csv`](mutedata.csv) | TCGA-LUAD mutation data | `1.get_data.R`: export KRAS/EGFR mutation classifications for duplicate-sample review. The first CSV column stores the original mutation-record barcode. |
| [`mutedata_modified.csv`](mutedata_modified.csv) | Manual curation of `mutedata.csv` | `1.get_data.R` reads this table to retain one mutation barcode per duplicated sample. This is a manually curated input, not a fully automatic output. |
| [`kras_sample.csv`](kras_sample.csv) | KRAS-mutant TCGA-LUAD expression and ImmuCellAI signatures | [`2_1.bulk_RNA_seq_analysis.R`](../Scripts/2_1.bulk_RNA_seq_analysis.R): compute ssGSEA scores from log2(TPM + 1), hierarchically cluster samples using the row-scaled score heatmap, split the sample tree with `cutree(k=2)`, and label the clusters cold/hot. Supplied table: 66 cold and 72 hot samples. |
| [`pde4d_tcga_30_samples.csv`](pde4d_tcga_30_samples.csv) | KRAS-mutant TCGA-LUAD PDE4D expression | Related grouping code is in [`4.PDE4D_analysis.R`](../Scripts/4.PDE4D_analysis.R): high = expression at or above the 70th percentile; low = expression at or below the 30th percentile; exclude intermediate samples. The supplied table has 42 high and 42 low samples, with sample IDs in the first CSV column. |
| [`pde4d_ccle_cutoff_tpmlogp1.csv`](pde4d_ccle_cutoff_tpmlogp1.csv) | KRAS-mutant LUAD cell lines from DepMap/CCLE **25Q3** | The LUAD/KRAS-mutant cell-line summary was assembled using mutation information. The PDE4D expression cutoff was set at 2.985 [log2(TPM + 1)], corresponding to a local minimum in the expression distribution. Expression > 2.985 is high; expression <= 2.985 is low. All supplied labels agree: **16 high and 12 low**; columns are `sample`, `expression`, and `group`. |

## Differential-expression tables

These CSVs contain full gene-level DESeq2 results, rather than only significant genes. Positive `log2FoldChange` means higher expression in the first group listed in the contrast column.

| File | Source / contrast | Generating script / method |
| --- | --- | --- |
| [`all_deseq2_results_mut_vs_wt.csv`](all_deseq2_results_mut_vs_wt.csv) | TCGA-LUAD: KRAS-mutant versus wild-type | `2_1.bulk_RNA_seq_analysis.R`: protein-coding counts, design `~ type`, with wild-type as the reference. |
| [`all_deseq2_results_cold_vs_hot.csv`](all_deseq2_results_cold_vs_hot.csv) | KRAS-mutant TCGA-LUAD: immune-cold versus immune-hot | `2_1.bulk_RNA_seq_analysis.R`: protein-coding counts, design `~ subtype`, with hot as the reference; grouping is recorded in `kras_sample.csv`. |
| [`all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv`](all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv) | KRAS-mutant TCGA-LUAD: PDE4D-high versus PDE4D-low | `4.PDE4D_analysis.R`: upper/lower expression-quantile groups; retain genes with total count > 10; design `~ group`; explicit contrast `high` versus `low`. |
| [`all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv`](all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv) | DepMap/CCLE KRAS-mutant LUAD cell lines: PDE4D-high versus PDE4D-low | `4.PDE4D_analysis.R`: density-cutoff groups and raw counts from default model entries; retain genes with total count > 10; design `~ group`; explicit contrast `high` versus `low`. |
| [`GSE183516_h358_ars1620_24h_vs_control_deseq2.csv`](GSE183516_h358_ars1620_24h_vs_control_deseq2.csv) | GEO GSE183516: H358, ARS-1620 at 24 h versus control | [`3.intersect_bulk_and_singlecell.R`](../Scripts/3.intersect_bulk_and_singlecell.R): select the six GEO sample columns named in the script, join the supplied GRCh38 annotation, retain genes with total count > 10, and fit DESeq2 with design `~ group` and contrast `24h` versus `control`. |
| [`GSE275534_h358_adag_72h_vs_control_deseq2.csv`](GSE275534_h358_adag_72h_vs_control_deseq2.csv) | GEO GSE275534: H358, adagrasib at 72 h versus control | `3.intersect_bulk_and_singlecell.R`: combine featureCounts files, select the eight columns used by the script, retain genes with total count > 10, and fit DESeq2 with design `~ group` and contrast `72h` versus `control`; map Ensembl gene IDs to HGNC symbols using BioMart. |

Common columns are:

- `gene_name`, `genes`, `gene`, or `hgnc_symbol`: gene identifier/symbol, as appropriate to the source table;
- `baseMean`: mean normalized count across samples in the DESeq2 analysis;
- `log2FoldChange`: estimated log2 fold change for the stated contrast;
- `lfcSE`: standard error of the log2 fold change;
- `stat`: Wald statistic;
- `pvalue`: unadjusted P value;
- `padj`: adjusted P value.

Significance filtering is applied downstream for specific analyses. For example, script `3` integrates bulk genes using `padj < 0.05` and absolute log2 fold change > log2(1.5); script `4` uses the same thresholds for the PDE4D KEGG analysis. The full result tables retain genes outside these thresholds.

## Integration, single-cell, and enrichment outputs

| File | Source | Generating script / method |
| --- | --- | --- |
| [`immune_genes_basical_info.csv`](immune_genes_basical_info.csv) | TCGA cold/hot DESeq2 results and public single-cell expression/immune-abundance correlations | `3.intersect_bulk_and_singlecell.R`: intersect cold-upregulated genes with negatively correlated single-cell genes, and hot-upregulated genes with positively correlated genes; single-cell selection uses correlation sign and unadjusted `p < 0.05`. Merge bulk/correlation statistics and BioMart gene descriptions; label immune association as negative or positive. Supplied table: 65 genes, 53 negativev(immune-suppressive) and 12 positive (immune-promoting). |
| [`overlap_65_genes_of_bulk_and_sc.csv`](overlap_65_genes_of_bulk_and_sc.csv) | The integrated CIIR gene table above | `3.intersect_bulk_and_singlecell.R`: add `nvalid_sample`, the count of nonmissing sample-level malignant-expression values, to the integrated table. |
| [`sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv`](sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv) | GSE136246, GSE123902, GSE146100, GSE148466, and GSE149655 | [`2_2.scRNA_seq_analysis.R`](../Scripts/2_2.scRNA_seq_analysis.R): summarize scaled malignant-cell expression and immune-cell fractions per sample; compute across-sample Pearson correlations; calculate BH-adjusted P values in `p.adj`. Supplied table: 59,427 rows; comma-separated sample-level vectors have 22 entries. This file is present locally (36.13 MiB); no current ignore rule excludes it. |
| [`umap_df_samples.csv`](umap_df_samples.csv) | The integrated public single-cell datasets above | `2_2.scRNA_seq_analysis.R`: export UMAP coordinates from the final Seurat object and attach cell identifiers, sample labels, and cell-type annotations. Columns: `UMAP_1`, `UMAP_2`, `cell`, `sample2`, `celltype`. Supplied table: 36,662 cells across 22 sample labels, including 31,708 `immun`, 3,226 `malignant`, and 1,728 `stromal` cells. |
| [`DEGs_enrichkegg_res_tcga_PDE4D.xlsx`](DEGs_enrichkegg_res_tcga_PDE4D.xlsx) | TCGA PDE4D-high versus low DESeq2 results | `4.PDE4D_analysis.R`: select up/down genes using `padj < 0.05` and absolute log2 fold change > log2(1.5), run KEGG enrichment through the script's `get_KEGG()` helper, and export up/down result tables. |

## External exports

These tables are consumed by the R scripts; their calculation/export is not implemented as a complete workflow in the current scripts. The mutation-query and uploaded-matrix branches are described separately below.

| File | Source / use | Input and export provenance |
| --- | --- | --- |
| [`TIMER2_KRAS_mutation_table.csv`](TIMER2_KRAS_mutation_table.csv) | TIMER2.0 mutation-associated immune-infiltration comparison; visualized in `2_1.bulk_RNA_seq_analysis.R` | We used TIMER2's built-in KRAS mutation query for this branch, without uploading a custom expression matrix. Columns include `cancer`, `infiltrates`, `log2FC`, `p`, and `adj.p`. |
| [`TIMER2_estimation_matrix_KRAStpm.csv`](TIMER2_estimation_matrix_KRAStpm.csv) | TIMER2.0 immune-infiltration estimates used by `4.PDE4D_analysis.R` | We uploaded a **TPM matrix for 138 KRAS-mutant TCGA-LUAD tumours**. The script prepares this input as `kras_exp.csv`. |

## Clinical output

[`integrated_data_49samples.csv`](integrated_data_49samples.csv) is a patient-level output of [`5.IHC_EFS_survival_curv.R`](../Scripts/5.IHC_EFS_survival_curv.R), merging institutional EFS and IHC inputs by `ID`.

It contains 49 records: 30 PDE4D-high and 19 PDE4D-low. The EFS, event, PDE4D H-score, CD8 percentage, group labels, and ID set match the [source-data workbook](<../Tables/Supplementary Data 8. Source data underlying the main and Supplementary figures.xlsx>), sheet `Fig6b-d`.

The manuscript uses these **49 KRAS-mutant LUAD specimens** for independent IHC support of the inverse association between PDE4D and **CD8⁺ T-cell abundance**, alongside the EFS analysis.

`EFS` is plotted in months; `Event` is the 0/1 event indicator passed to `Surv()`. `PDE4D` is the H-score, `CD8_%` is the CD8 percentage, and `group` is high for H-score > 50 and low otherwise. The original script uses `CD8` as a column name after raw-input preparation; when loading this processed CSV directly, adapt that name accordingly.
## Data provenance and grouping notes

The [main README](../README.md#data-sources-and-required-inputs) records dataset accessions, sources, and dataset releases/download information.

The author reports that most analyses used R 4.3.1; some scripts may require older R/package versions, including the `pheatmap` environment used for hot/cold grouping. See [Reproducibility status](../README.md#reproducibility-status) for the remaining filename and execution-order issues in the current scripts.
