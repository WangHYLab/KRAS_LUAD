# Processed results and sample metadata

This directory contains processed tables used by the analyses described in the [main README](../README.md). Raw-data access and dataset-version information are documented there.

**Terminology:** CIIR genes means **cancer-cell-intrinsic immune-related (CIIR) genes in KRAS-mutant LUAD**. The integrated gene tables below contain the CIIR gene analyses and associated statistics.

The descriptions below are based on the current scripts and the supplied table headers. The tables have not been regenerated or compared numerically with a new analysis run. Where an exact export statement is absent, that is stated explicitly rather than assuming a verified generation history.

## Sample lists and grouping tables

| File | Source | Generating script / method |
| --- | --- | --- |
| [`sk.csv`](sk.csv) | TCGA-LUAD expression and mutation data | [`1.get_data.R`](../Scripts/1.get_data.R): match expression and mutation barcodes using their first 16 characters, apply the recorded duplicate-barcode selection, retain primary-tumour barcodes ending in `01A`, and assign KRAS status. Supplied table: 502 samples, with `KRAS=1` for 138 mutant and `KRAS=0` for 364 wild-type samples. |
| [`mutedata.csv`](mutedata.csv) | TCGA-LUAD mutation data | `1.get_data.R`: export KRAS/EGFR mutation classifications for duplicate-sample review. The first CSV column stores the original mutation-record barcode. |
| [`mutedata_modified.csv`](mutedata_modified.csv) | Manual curation of `mutedata.csv` | `1.get_data.R` reads this table to retain one mutation barcode per duplicated sample. This is a manually curated input, not a fully automatic output. [TO COMPLETE: selection rationale]. |
| [`kras_sample.csv`](kras_sample.csv) | KRAS-mutant TCGA-LUAD expression and ImmuCellAI signatures | [`2_1.bulk_RNA_seq_analysis.R`](../Scripts/2_1.bulk_RNA_seq_analysis.R): compute ssGSEA scores from log2(TPM + 1), hierarchically cluster samples using the row-scaled score heatmap, split the sample tree with `cutree(k=2)`, and label the clusters cold/hot. Supplied table: 66 cold and 72 hot samples. |
| [`pde4d_tcga_30_samples.csv`](pde4d_tcga_30_samples.csv) | KRAS-mutant TCGA-LUAD PDE4D expression | Related grouping code is in [`4.PDE4D_analysis.R`](../Scripts/4.PDE4D_analysis.R): high = expression at or above the 70th percentile; low = expression at or below the 30th percentile; exclude intermediate samples. The supplied table has 42 high and 42 low samples, with sample IDs in the first CSV column. The current script does not explicitly export this filename; confirm its correspondence to the final selected-sample object. |
| [`pde4d_ccle_cutoff_tpmlogp1.csv`](pde4d_ccle_cutoff_tpmlogp1.csv) | KRAS-mutant LUAD cell lines from DepMap/CCLE | Related code in `4.PDE4D_analysis.R` applies a density-valley cutoff to PDE4D log2(TPM + 1), with expression above the cutoff labelled high and the remainder low. Supplied table: 16 high and 12 low cell lines; columns are `sample`, `expression`, and `group`. The script exports `samples_PDE4D_group_cutoff_valley_ccle.csv`, rather than this filename. [TO COMPLETE: verify the saved table's exact cutoff and export history]. |

The names of the PDE4D group files should not be used to infer sample counts. In particular, `pde4d_tcga_30_samples.csv` contains 84 samples; the related script uses upper/lower 30% expression groups.

## Differential-expression tables

These CSVs contain full gene-level DESeq2 results, rather than only significant genes. Positive `log2FoldChange` means higher expression in the first group listed in the contrast column.

| File | Source / contrast | Generating script / method |
| --- | --- | --- |
| [`all_deseq2_results_mut_vs_wt.csv`](all_deseq2_results_mut_vs_wt.csv) | TCGA-LUAD: KRAS-mutant versus wild-type | `2_1.bulk_RNA_seq_analysis.R`: protein-coding counts, design `~ type`, with wild-type as the reference. The script also reloads `dds2.RData` before exporting; the supplied result's checkpoint version should be recorded. |
| [`all_deseq2_results_cold_vs_hot.csv`](all_deseq2_results_cold_vs_hot.csv) | KRAS-mutant TCGA-LUAD: immune-cold versus immune-hot | `2_1.bulk_RNA_seq_analysis.R`: protein-coding counts, design `~ subtype`, with hot as the reference; grouping is recorded in `kras_sample.csv`. |
| [`all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv`](all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv) | KRAS-mutant TCGA-LUAD: PDE4D-high versus PDE4D-low | `4.PDE4D_analysis.R`: upper/lower expression-quantile groups; retain genes with total count > 10; design `~ group`; explicit contrast `high` versus `low`. |
| [`all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv`](all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv) | DepMap/CCLE KRAS-mutant LUAD cell lines: PDE4D-high versus PDE4D-low | `4.PDE4D_analysis.R`: density-cutoff groups and raw counts from default model entries; retain genes with total count > 10; design `~ group`; explicit contrast `high` versus `low`. |
| [`GSE183516_h358_ars1620_24h_vs_control_deseq2.csv`](GSE183516_h358_ars1620_24h_vs_control_deseq2.csv) | GEO GSE183516: H358, ARS-1620 at 24 h versus control | [`3.intersect_bulk_and_singlecell.R`](../Scripts/3.intersect_bulk_and_singlecell.R): select the six GEO sample columns named in the script, join the supplied GRCh38 annotation, retain genes with total count > 10, and fit DESeq2 with design `~ group` and contrast `24h` versus `control`. |
| [`GSE275534_h358_adag_72h_vs_control_deseq2.csv`](GSE275534_h358_adag_72h_vs_control_deseq2.csv) | GEO GSE275534: H358, adagrasib at 72 h versus control | `3.intersect_bulk_and_singlecell.R`: combine featureCounts files, select the eight columns used by the script, retain genes with total count > 10, and fit DESeq2 with design `~ group` and contrast `72h` versus `control`; map Ensembl gene IDs to HGNC symbols using BioMart. The current write statement uses a filename without `h358`; align the export/read names when reproducing this table. |

Common columns are:

- `gene_name`, `genes`, `gene`, or `hgnc_symbol`: gene identifier/symbol, as appropriate to the source table;
- `baseMean`: mean normalized count across samples in the DESeq2 analysis;
- `log2FoldChange`: estimated log2 fold change for the stated contrast;
- `lfcSE`: standard error of the log2 fold change;
- `stat`: Wald statistic;
- `pvalue`: unadjusted P value;
- `padj`: adjusted P value. Missing values may occur and should not be replaced with zero.

Significance filtering is applied downstream for specific analyses. For example, script `3` integrates bulk genes using `padj < 0.05` and absolute log2 fold change > log2(1.5); script `4` uses the same thresholds for the PDE4D KEGG analysis. The full result tables retain genes outside these thresholds.

## Integration, single-cell, and enrichment outputs

| File | Source | Generating script / method |
| --- | --- | --- |
| [`immune_genes_basical_info.csv`](immune_genes_basical_info.csv) | TCGA cold/hot DESeq2 results and public single-cell expression/immune-abundance correlations | `3.intersect_bulk_and_singlecell.R`: intersect cold-upregulated genes with negatively correlated single-cell genes, and cold-downregulated genes with positively correlated genes; single-cell selection uses correlation sign and unadjusted `p < 0.05`. Merge bulk/correlation statistics and BioMart gene descriptions; label immune association as negative or positive. |
| [`overlap_65_genes_of_bulk_and_sc.csv`](overlap_65_genes_of_bulk_and_sc.csv) | The integrated CIIR gene table above | `3.intersect_bulk_and_singlecell.R`: add `nvalid_sample`, the count of nonmissing sample-level malignant-expression values, to the integrated table. The number in the filename is a label for the saved analysis and is not a general filtering rule. |
| `sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv` | GSE136246, GSE123902, GSE146100, GSE148466, and GSE149655 | [`2_2.scRNA_seq_analysis.R`](../Scripts/2_2.scRNA_seq_analysis.R): summarize scaled malignant-cell expression and immune-cell fractions per sample; compute across-sample Pearson correlations; calculate BH-adjusted P values in `p.adj`. Includes comma-separated sample-level values. This large table is currently excluded by `.gitignore`; supply it through an archive or regenerate it. |
| [`umap_df_samples.csv`](umap_df_samples.csv) | The integrated public single-cell datasets above | `2_2.scRNA_seq_analysis.R`: export UMAP coordinates from the final Seurat object and attach cell identifiers, sample labels, and cell-type annotations. Columns: `UMAP_1`, `UMAP_2`, `cell`, `sample2`, `celltype`. |
| [`DEGs_enrichkegg_res_tcga_PDE4D.xlsx`](DEGs_enrichkegg_res_tcga_PDE4D.xlsx) | TCGA PDE4D-high versus low DESeq2 results | `4.PDE4D_analysis.R`: select up/down genes using `padj < 0.05` and absolute log2 fold change > log2(1.5), run KEGG enrichment through the script's `get_KEGG()` helper, and export up/down result tables. |

## External exports

These tables are consumed by the R scripts; their calculation/export is not implemented as a complete workflow in the current scripts.

| File | Source / use | Information to record |
| --- | --- | --- |
| [`TIMER2_KRAS_mutation_table.csv`](TIMER2_KRAS_mutation_table.csv) | TIMER2.0 mutation-associated immune-infiltration comparison; visualized in `2_1.bulk_RNA_seq_analysis.R` | [TO COMPLETE: export date, cancer/mutation query, tool settings, and resource version]. Columns include `cancer`, `infiltrates`, `log2FC`, `p`, and `adj.p`. |
| [`TIMER2_estimation_matrix_KRAStpm.csv`](TIMER2_estimation_matrix_KRAStpm.csv) | TIMER2.0 immune-infiltration estimates used by `4.PDE4D_analysis.R` | The script prepares the KRAS-mutant TCGA TPM matrix as `kras_exp.csv`. [TO COMPLETE: confirm the uploaded matrix and record export date, selected methods, settings, and expression units]. The supplied table has immune-estimate rows and 138 sample columns. |

## Clinical output

[`integrated_data_49samples.csv`](integrated_data_49samples.csv) is a patient-level output of [`5.IHC_EFS_survival_curv.R`](../Scripts/5.IHC_EFS_survival_curv.R), merging institutional EFS and IHC inputs by `ID`. The author has confirmed that the processed IHC/EFS data can be publicly shared, and this file is included in the repository.

It contains 49 records. The EFS, event, PDE4D H-score, CD8 percentage, group labels, and ID set match the [source-data workbook](<../Tables/Supplementary Data 9. Source data underlying the main and Supplementary figures.xlsx>), sheet `Fig6b-d`. This is an institutional clinical table, distinct from the public TCGA/GEO-derived results.

`EFS` is plotted in months; `Event` is the 0/1 event indicator passed to `Surv()`. `PDE4D` is the H-score, `CD8_%` is the CD8 percentage, and `group` is high for H-score > 50 and low otherwise. The original script uses `CD8` as a column name after raw-input preparation; when loading this processed CSV directly, adapt that name accordingly. Original raw clinical input CSVs are not bundled.

## Metadata to complete before release

Record the raw-data versions/download dates, the original R/package environment, the manually curated barcode-selection rationale, and the exact export history/cutoff of the saved PDE4D grouping tables. External TIMER2 exports need their own provenance details. See the main README for the manuscript citation, data archive, and contact information.

Most analyses used R 4.3.1; some scripts may require older R/package versions, including the `pheatmap` environment used for hot/cold grouping. The author recalls DepMap/CCLE release 25Q3, pending verification.
