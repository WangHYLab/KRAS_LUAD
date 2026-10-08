# Integrative analysis reveals immune heterogeneity and cancer-cell-intrinsic immune-related genes in KRAS-mutant lung adenocarcinoma

Analysis code accompanying a study of the tumour immune microenvironment and tumour-cell programs associated with immune heterogeneity in KRAS-mutant lung adenocarcinoma (LUAD). The study identifies immune-hot and immune-cold subtypes and prioritizes PDE4D as a KRAS-responsive cancer-cell-intrinsic immune-related (CIIR) candidate associated with reduced infiltration across multiple immune-cell populations, alterations in cytokine-related transcriptional programs, and lower tumour-cell HLA-I expression.

**Terminology:** CIIR refers to **cancer-cell-intrinsic immune-related genes in KRAS-mutant LUAD**. The integrated bulk/single-cell tables identify CIIR genes; KRAS responsiveness is assessed using KRAS-inhibitor perturbation data. PDE4D is described as a **KRAS-responsive CIIR candidate**. Immune-hot/immune-cold denotes tumour immune subtypes; PDE4D-high/PDE4D-low denotes the expression or IHC groups defined for each analysis. Gene associations are described as positive or negative with respect to immune-cell abundance, matching the saved `positive`/`negative` labels.

The analyses combine TCGA-LUAD bulk RNA-sequencing and mutation data, public single-cell RNA-sequencing datasets, KRAS-inhibitor perturbation datasets, DepMap/CCLE cell-line data, and tissue immunohistochemistry (IHC) and event-free survival (EFS) analyses.

The manuscript also reports enhanced tumour-cell suppression after **PDE4 inhibition** in whole-blood co-culture and increased **MHC-I components** after pharmacological and genetic perturbation. These experimental findings are described at the study level; the computational analyses and IHC/EFS workflow available here are mapped below.

**Publication status:** Manuscript under submission. Publication details and a DOI will be added when available.

**Release status:** Local research-code collection; not yet a verified reproducible release. The scripts require input preparation, path configuration, and the corrections described under [Reproducibility status](#reproducibility-status). The full workflow has not been verified in a clean environment.

**Documentation check (2026-10-08):** This README reflects the files currently present, the current script references, and counts checked directly from the supplied tables. It does not establish that the scripts regenerate those tables or the manuscript figures.

## Repository contents

```text
Scripts/       Six R analysis scripts
Tables/        Immune signatures, reference tables, and source-data workbook
Results/       Processed analysis tables and sample metadata
Rdata/         Saved R objects available locally
Figures/       Assembled manuscript figure PDFs
README.md      Study overview, inputs, and workflow
```

The current directory is not initialized as a Git repository, and neither `.gitignore` nor `LICENSE` is present. No Git tracking or exclusion status is asserted here. `Rdata/`, `Figures/`, the full single-cell correlation CSV, the processed IHC/EFS CSV, and the source-data workbook are all present locally. The contents of a future public release will depend on the files actually published.

## Analysis workflow and figure mapping

Figure mappings follow comments in the current scripts and should be checked against the final manuscript.

| Script | Purpose | Associated figures |
| --- | --- | --- |
| `Scripts/1.get_data.R` | Download TCGA-LUAD expression and mutation data; reconcile sample barcodes; prepare matched cohorts | Input preparation |
| `Scripts/2_1.bulk_RNA_seq_analysis.R` | KRAS-mutant versus wild-type differential expression; immune-signature GSEA; ssGSEA immune subtyping; mutation, survival, and functional comparisons | Figures 2–3; Supplementary Figure 1 |
| `Scripts/2_2.scRNA_seq_analysis.R` | Single-cell preprocessing, CopyKAT inference, integration and annotation; correlate malignant-cell expression with immune-cell abundance | Supplementary Figure 2; inputs to Figure 4 |
| `Scripts/3.intersect_bulk_and_singlecell.R` | Integrate bulk and single-cell associations to identify CIIR genes; annotate shared genes; assess KRAS responsiveness using KRAS-inhibitor perturbation data | Figure 4b–h |
| `Scripts/4.PDE4D_analysis.R` | PDE4D stratification; associations with infiltration across multiple immune-cell populations; cytokine-related transcriptional programs and KEGG enrichment; HLA-I, NLRC5 and CREB1 expression; inferred JAK-STAT pathway and STAT1 activity | Figure 5; Figure 6e–g; Supplementary Figures 3–5 and 7 |
| `Scripts/5.IHC_EFS_survival_curv.R` | PDE4D IHC scoring, associations with CD8⁺ T-cell abundance, and EFS analysis | Figure 6b–d |

The current script set does not include a dedicated workflow for the study schematic, conceptual model, CREB1 ChIP-seq panel, or every experimental/image panel. The repository describes the computational analyses present here; it does not claim to regenerate every manuscript panel automatically.

## Data sources and required inputs

| Source | Role | Access / preparation |
| --- | --- | --- |
| TCGA-LUAD | Bulk expression and somatic mutations | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/); downloaded in **2024**, as confirmed by the author; queried with `TCGAbiolinks` in `1.get_data.R` |
| GSE136246, GSE123902, GSE146100, GSE148466, GSE149655 | Single-cell analysis | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/); **22 KRAS-mutant samples** selected using the sample information provided in the GEO datasets; accession-specific input files are referenced in `2_2.scRNA_seq_analysis.R` |
| GSE183516 | H358 ARS-1620, 24 h versus control, as selected in the script | [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183516); count matrix and gene annotation table |
| GSE275534 | H358 adagrasib, 72 h versus control, as selected in the script | [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275534); featureCounts files |
| DepMap / CCLE | KRAS-mutant LUAD cell-line expression and PDE4D grouping | [DepMap downloads](https://depmap.org/portal/download/), release **25Q3**, confirmed by the author; the LUAD/KRAS-mutant cell-line summary was assembled using the mutation information |
| TIMER2.0 | Immune-infiltration estimates and KRAS-associated comparisons | [TIMER2.0](http://timer.cistrome.org/); analyses were performed during **2025–2026**. The estimation branch used an uploaded TPM matrix for **138 KRAS-mutant tumours**; mutation-associated comparisons used the built-in KRAS mutation query without uploading a custom matrix |
| CIBERSORT signatures | Immune gene-set enrichment | CIBERSORT **LM22 immune cell signatures**; the marker table used by the scripts is `Tables/CIBERSORT_cell_types_markers.csv` |
| ImmuCellAI signatures | Immune gene-set enrichment and ssGSEA subtyping | Immune-cell signatures provided by the **ImmuCellAI tool**; the supplied table is `Tables/ImmuCellAI_24_immune_cell_types_signature_genes.csv` |
| TCGA survival table | Survival comparison of immune subtypes | Downloaded from **UCSC Xena**; supplied as `Tables/survival_LUAD_survival.txt`. Overall-survival time (`OS.time`) is in **days**; `OS` is the event indicator |
| Institutional IHC / clinical cohort | PDE4D–CD8⁺ T-cell associations and EFS analyses | Local source data are provided in the source-data workbook, sheet `Fig6b-d`, and `Results/integrated_data_49samples.csv`; the original input CSVs read by script `5` are not included |

**Dataset provenance:** The author confirmed the TCGA download year (**2024**), DepMap/CCLE release (**25Q3**), UCSC Xena survival-data source, and signature sources above. TIMER2 was used in **2025–2026**; exact export dates for the individual tables are not recorded. Public databases can change; accession IDs alone do not freeze the exact analysis inputs.

For the DepMap/CCLE analysis, a PDE4D expression cutoff of **2.985** on the **log2(TPM + 1)** scale, corresponding to a local minimum of the expression distribution, classified **16 cell lines as PDE4D-high** and **12 as PDE4D-low**. The script uses expression > cutoff for high and expression <= cutoff for low; all 28 labels in the supplied grouping table agree with this rule.

The current TCGA sample metadata contain 502 tumours: 138 KRAS-mutant and 364 KRAS-wild-type. The KRAS-mutant subtype table contains 66 immune-cold and 72 immune-hot tumours. These counts describe the supplied metadata, rather than a guarantee that a new database query will return identical cohorts.

The TCGA preprocessing script uses a manually curated mutation-barcode selection (`Results/mutedata_modified.csv`) and a fixed positional selection of duplicate expression barcodes. The supplied table and current code retain these selections; the original selection rationale is not recorded.

### Processed outputs

See [Results/README.md](Results/README.md) for a file-by-file description of source datasets, generating scripts, grouping rules, and differential-expression contrasts. Paths in the scripts are intentionally anonymized; users should configure their own local directories.

Selected local outputs include:

- `Results/all_deseq2_results_mut_vs_wt.csv` and `Results/all_deseq2_results_cold_vs_hot.csv`;
- `Results/sk.csv`, `Results/kras_sample.csv`, and `Results/mutedata_modified.csv`;
- `Results/overlap_65_genes_of_bulk_and_sc.csv` and `Results/immune_genes_basical_info.csv`: each contains 65 genes, labelled negative (53) or positive (12) for immune association;
- `Results/sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv`: 59,427 rows with sample-level vectors of length 22;
- `Results/umap_df_samples.csv`: 36,662 cells from 22 sample labels;
- KRAS-inhibitor differential-expression tables for GSE183516 and GSE275534;
- PDE4D group tables, differential-expression results, and TIMER2 exported estimates.

`Rdata/` contains TCGA expression/mutation matrices, the cold/hot DESeq2 object, PDE4D DESeq2 objects, and several enrichment/activity objects. The full single-cell correlation CSV is also present (36.13 MiB). The local `dds_PDE4D_tcga.RData` is 265.46 MiB. These file sizes need to be considered when preparing a code/data release; local presence does not establish public download availability.

Some required inputs are absent, including `dds2.RData`, the original public-data downloads, single-cell Seurat/RDS objects, and the original clinical CSVs. Some other files absent from this directory can be generated by earlier script sections, as described below. Before a reproducible release, document how each required input is regenerated or provide it through an archive with filenames and checksums.

**Processed-data archive:** No persistent archive URL/DOI or checksum manifest is currently documented.

### Clinical data

The IHC/EFS script reads a clinical EFS CSV and `PDE4D_CD8.csv`, merges them by `ID`, and uses a PDE4D H-score cutoff of 50: scores > 50 are labelled high, and scores <= 50 are labelled low.

The manuscript uses this independent IHC cohort of **49 KRAS-mutant LUAD specimens** to support the inverse association between PDE4D and **CD8⁺ T-cell abundance**.

The [source-data workbook](<Tables/Supplementary Data 9. Source data underlying the main and Supplementary figures.xlsx>), sheet `Fig6b-d`, contains 49 records with `ID`, `EFS`, `Event`, `PDE4D_Percentage`, `PDE4D_intensity`, `PDE4D_Hscore`, `CD8_%`, and `PDE4D_group`. Its IDs and EFS, event, H-score, CD8, and group values match [the processed IHC/EFS CSV](Results/integrated_data_49samples.csv).

The supplied processed table contains 30 PDE4D-high and 19 PDE4D-low records under the H-score cutoff above.

**Clinical data availability:** The processed IHC/EFS CSV and the source-data workbook are included for public sharing, as confirmed by the author. The original clinical input CSVs are not bundled; the supplied processed table contains the measurements and group labels underlying Figure 6b–d.

## Software environment

The author reports that most analyses used **R 4.3.1**. **Some scripts may require older versions of R or their dependent packages**, particularly the immune-hot/immune-cold clustering, which used an older `pheatmap` version.

Package references observed in the scripts include:

| Area | Packages |
| --- | --- |
| Data acquisition and annotation | `TCGAbiolinks`, `scRNAseq`, `biomaRt`, `org.Hs.eg.db` |
| Differential expression, enrichment, mutation analysis | `DESeq2`, `clusterProfiler`, `GSVA`, `GseaVis`, `maftools` |
| Single-cell analysis | `Seurat`, `copykat`, `future` |
| Signalling and transcription-factor activity | `progeny`, `dorothea`, `viper` |
| Data processing | `tidyverse`, `dplyr`, `tidyr`, `tibble`, `readr`, `purrr`, `stringr`, `reshape2`, `openxlsx`, `mgcv` |
| Plotting and survival | `ggplot2`, `ggpubr`, `ggrepel`, `pheatmap`, `ComplexHeatmap`, `patchwork`, `gridExtra`, `viridis`, `survival`, `survminer` |

Use the original analysis environment where possible. Some scripts use older function interfaces, including the legacy `GSVA::gsva()` call, and saved objects can depend on package versions. Install dependencies using their official distribution sources and use compatible versions where needed.

The reported analysis environment has not been independently rerun during repository preparation.

## Preparing and running the analyses

1. Decide whether to recompute from raw public data or reuse the supplied processed tables/local R objects. Obtain the missing inputs for the chosen analysis sections; no processed-data archive link is currently provided.
2. Replace the intentional placeholder paths (`path/to/...`, `/xx/...`, and similar examples) with your own local paths. The author's private filesystem paths are not needed. Configure `WORK_DIR`, `setwd()`, and the input/output paths in each script. Use forward slashes in R paths and a trailing slash where `paste0(WORK_DIR, ...)` expects one. A single root-directory change does not currently update all script paths.
3. Resolve the filename, execution-order, and dependency issues listed below. Create the output directories referenced by the scripts. Bare filenames are resolved relative to `setwd()`, whereas other sections use `WORK_DIR`; configuring one does not configure the other.
4. For raw-data recomputation, run `1.get_data.R` for TCGA input preparation. Complete its duplicate-barcode selections before continuing. Its relative output paths must be aligned with the `Results/` and `Rdata/` locations used downstream.
5. Run the bulk (`2_1`) and single-cell (`2_2`) branches, then the integration and perturbation analysis (`3`). Run the PDE4D analysis (`4`) once its TCGA, signature, TIMER2, and cell-line inputs are available.
6. Run the IHC/EFS branch (`5`) separately. Its initial preparation section expects the original clinical input CSVs; use the supplied processed IHC/EFS table when reproducing downstream analyses, adapting the data-loading and column-name steps as needed.

After resolving the issues below, run the relevant sections in R/RStudio. The scripts mix computation, saved-object reloads, and plotting, and are not yet packaged as a single command-line pipeline. Some plots are displayed as R objects; assembled manuscript PDFs are separate artifacts.

## Reproducibility status

The following issues remain in the current files. Placeholder paths are intentional and should be configured locally.

| Location | Current issue / action needed |
| --- | --- |
| `1.get_data.R` and downstream scripts | The preparation script writes files to its working directory; downstream sections mix bare filenames with `WORK_DIR/Results/` and `WORK_DIR/Rdata/`. Align these paths before running. |
| `2_1.bulk_RNA_seq_analysis.R`, line 61 | `load('dds2.RData')` unconditionally replaces the newly fitted KRAS-mutant/wild-type object, but that checkpoint is absent and its save statement is commented out. Establish the intended fit/checkpoint before exporting. |
| `3.intersect_bulk_and_singlecell.R`, lines 94 and 139 | The bulk heatmap uses `data5` before the single-cell heatmap section creates it. Correct the execution order and confirm the intended gene order. |
| `3.intersect_bulk_and_singlecell.R`, line 99 | The single-cell heatmap reads a concatenated correlation filename that is absent. Confirm the intended table; the standard `sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv` is supplied. |
| `3.intersect_bulk_and_singlecell.R`, lines 323/327 | GSE275534 results are written without `h358` in the name and read with it. Align the export/read names with the supplied result. |
| `3.intersect_bulk_and_singlecell.R`, lines 331/397 | The script reads `fig4.overlap_65_genes_of_bulk_and_sc.csv`; the supplied and earlier exported file is `overlap_65_genes_of_bulk_and_sc.csv`. Confirm the analysis version before aligning names. |
| `4.PDE4D_analysis.R`, line 1590 | The CIIR annotation table is read from `Tables/immune_genes_basical_info.csv`, while the supplied table is in `Results/`. Align the directory reference. |
| TCGA selection, clustering, and plotting sections | Expression-column ranges, duplicate-barcode positions, sample counts, cluster labels, and some plotting selections are fixed. Revalidate these choices when inputs or package versions change. |
| Single-cell and clinical branches | Raw inputs and single-cell RDS objects are absent. Regenerate the single-cell objects from the public inputs; use the supplied clinical table for downstream plots with appropriate column-name/loading changes. |

The obsolete filenames `all_deseq2_results_cold_vs_hot_pc.csv`, `ImmuCellAI_immucellSignature.csv`, and `CIBERSORT_22_immune_cell_types_signature_genes.csv` are no longer referenced by the current scripts. Their names now match the supplied tables, although their directory paths still need configuration.

`Rdata/dds2_ch.RData` is present and is also generated by the bulk script. The CCLE grouping file `samples_PDE4D_group_cutoff_valley_ccle.csv` and `timer2_spearman_cor_plot_df.Rdata` are absent locally but have active generation statements earlier in script `4`; they are not separate missing external inputs when those sections run successfully. The supplied `pde4d_ccle_cutoff_tpmlogp1.csv` uses a different filename and includes an expression column; its labels agree with the author-confirmed cutoff of 2.985. The script computes a density-valley cutoff from its inputs rather than hard-coding that reported value.

`pde4d_tcga_30_samples.csv` contains 84 samples (42 per group), and the TCGA grouping code uses the upper and lower 30% of PDE4D expression. The current script does not explicitly export that grouping filename, so correspondence to the final saved analysis should still be checked.

This documentation update checked file presence, current script references, table counts, and agreement between the clinical CSV and workbook. No R analysis was rerun. Regeneration of results and figures, checkpoint contents, and compatibility with a clean R environment remain unverified.

Dataset sources and release/download information belong in the data-source table above. File-specific sample selection, grouping cutoffs, manual curation criteria, and external-export settings belong in [Results/README.md](Results/README.md#data-provenance-and-grouping-notes).

## Citation and contact

Please cite the associated manuscript when publicly available, together with the original datasets/methods when using these analyses.

**Manuscript:** *Integrative analysis reveals immune heterogeneity and cancer-cell-intrinsic immune-related genes in KRAS-mutant lung adenocarcinoma*. Under submission; author and publication details will be added later.

**Code release DOI/version:** Not yet documented.

**Corresponding author / contact:** To be added later.

## License

The author previously selected the **MIT License** for the analysis code. The current directory does not contain a `LICENSE` file; add the intended license text and copyright-holder name before publishing the code release.

The code license does not automatically apply to third-party datasets, signatures, manuscript figures, or source data; these materials retain their applicable terms.
