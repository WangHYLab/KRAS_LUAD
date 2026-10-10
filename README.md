# Integrative analysis reveals cancer-cell-intrinsic genes associated with immune heterogeneity in KRAS-mutant lung adenocarcinoma

Analysis code accompanying a study of the tumour immune microenvironment and tumour-cell programs associated with immune heterogeneity in KRAS-mutant lung adenocarcinoma (LUAD). The study identifies immune-hot and immune-cold subtypes and prioritizes PDE4D as a KRAS-responsive cancer-cell-intrinsic immune-related (CIIR) candidate associated with reduced infiltration across multiple immune-cell populations, alterations in cytokine-related transcriptional programs, and reduced HLA-I-related expression.

The analyses combine TCGA-LUAD bulk RNA-sequencing and mutation data, public single-cell RNA-sequencing datasets, KRAS-inhibitor perturbation datasets, DepMap cell-line data, and tissue immunohistochemistry (IHC) data.

**Publication status:** Manuscript under submission. Publication details and a DOI will be added when available.

## Repository contents

```text
Scripts/       Six R analysis scripts
Tables/        Immune signatures, reference tables, and source-data workbook
Results/       Processed analysis tables and sample metadata
README.md      Study overview, inputs, and workflow
```

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

## Data sources and required inputs

| Source | Role | Access / preparation |
| --- | --- | --- |
| TCGA-LUAD | Bulk expression and somatic mutations | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/); queried with `TCGAbiolinks` in `1.get_data.R` |
| GSE136246, GSE123902, GSE146100, GSE148466, GSE149655 | Single-cell analysis | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/); **22 KRAS-mutant samples** selected using the sample information provided in the GEO datasets; accession-specific input files are referenced in `2_2.scRNA_seq_analysis.R` |
| GSE183516 | H358 ARS-1620, 24 h versus control; 3 control samples and 3 treated samples, as selected for this analysis | [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183516); count matrix and gene annotation table |
| GSE275534 | H358 MRTX849 (adagrasib), 72 h versus control; 4 control samples and 4 treated samples, as selected for this analysis | [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275534); featureCounts files |
| DepMap / CCLE | KRAS-mutant LUAD cell-line expression and PDE4D grouping | [DepMap downloads](https://depmap.org/portal/download/), release **25Q3**; the LUAD/KRAS-mutant cell-line summary was assembled using the mutation information |
| TIMER2.0 | Immune-infiltration estimates and KRAS-associated comparisons | [TIMER2.0](http://timer.cistrome.org/); analyses were performed during **2025–2026**. Immune-infiltration estimation was performed using an uploaded TPM matrix for 138 KRAS-mutant tumours, whereas KRAS mutation-associated comparisons used the built-in mutation module without uploading a custom expression matrix |
| CIBERSORT signatures | Immune gene-set enrichment | CIBERSORT **LM22 immune cell signatures**; the LM22-derived signature table used by the scripts is `Tables/CIBERSORT_cell_types_markers.csv` |
| ImmuCellAI signatures | Immune gene-set enrichment and ssGSEA subtyping | Immune-cell signature genes derived from ImmuCellAI; the supplied table is `Tables/ImmuCellAI_24_immune_cell_types_signature_genes.csv` |
| TCGA survival table | Survival comparison of immune subtypes | Downloaded from **UCSC Xena**; supplied as `Tables/survival_LUAD_survival.txt`. Overall-survival time (`OS.time`) is in **days**; `OS` is the event indicator |
| Public CREB1 ChIP-seq | CREB1 binding at HLA-A, HLA-B and HLA-C loci; Supplementary Figure 6 | Cistrome DB ID **37305** (as identified by the author); [GEO GSM1010720](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1010720); raw sequencing data are linked under SRA **SRX190176**. Existing processed files were visualized with IGV; no new ChIP-seq data were generated |
| Institutional IHC / clinical cohort | PDE4D–CD8⁺ T-cell associations and EFS analyses | Local source data are provided in the source-data workbook, sheet `Fig6b-d`, and `Results/integrated_data_49samples.csv`; the original input CSVs read by script `5` are not included |

For the DepMap cancer cell lines analysis, a PDE4D expression cutoff of **2.985** on the **log2(TPM + 1)** scale, corresponding to a local minimum of the expression distribution, classified **16 cell lines as PDE4D-high** and **12 as PDE4D-low**. The script uses expression > cutoff for high and expression <= cutoff for low; all 28 labels in the supplied grouping table agree with this rule.

The current TCGA sample metadata contain 502 tumours: 138 KRAS-mutant and 364 KRAS-wild-type. The KRAS-mutant subtype table contains 66 immune-cold and 72 immune-hot tumours. These counts reflect the supplied analysis metadata and may not be reproduced exactly by a new query of the source databases.

The TCGA preprocessing script uses a manually curated mutation-barcode selection (Results/mutedata_modified.csv) and a fixed index-based selection among duplicate expression barcodes.

### Processed outputs

See [Results/README.md](Results/README.md) for a file-by-file description of processed datasets and result tables, generating scripts, grouping rules, and differential-expression contrasts.

### Clinical data

The IHC/EFS script reads a clinical EFS CSV and `PDE4D_CD8.csv`, merges them by `ID`, and uses a PDE4D H-score cutoff of 50: scores > 50 are labelled high, and scores <= 50 are labelled low.

The manuscript uses this independent IHC cohort of **49 KRAS-mutant LUAD specimens** to support the inverse association between PDE4D and **CD8⁺ T-cell abundance**.

## Software environment

The author reports that most analyses used **R 4.3.1**. Analyses were performed at different times; the versions listed below should not be interpreted as a single fully tested installation environment. Version entries are author-reported rather than reconstructed from a complete archived `sessionInfo()` record.

| Software/package | Reported version | Scope / version detail |
| --- | --- | --- |
| R | 4.3.1 | Reported for most analyses |
| DESeq2 | 1.42.0 | Differential-expression analysis |
| clusterProfiler | 4.10.0 | GSEA and GO/KEGG enrichment |
| org.Hs.eg.db | 3.18.0 | Human gene annotation |
| GSVA | **1.48 series** | Original ssGSEA analysis; the author confirmed the release series but the exact patch version could not be recovered |
| Seurat | 4.0.1 | Single-cell analysis |
| SingleR | 1.4.1 | Reported for single-cell annotation; confirm the corresponding use in the final script set |
| progeny | 1.24.0 | Pathway activity estimation |
| dorothea | 1.14.1 | Transcription-factor regulons |
| viper | 1.36.0 | Transcription-factor activity estimation |
| survival | 3.5.7 | Survival analysis |
| survminer | 0.4.9 | Survival visualization |
| pheatmap | Exact version not recovered | **1.0.12 is a candidate based on the recalled analysis period and compatible code; it has not been confirmed as the version actually used** |
| IGV | 2.19.2 | Visualization of public CREB1 ChIP-seq data |
| TIMER2.0 | Web application | Immune-infiltration analyses; access period is reported in the data-source table |
| ImageJ | Version to be confirmed | Immunoblot intensity quantification |

Other package references observed in the scripts include `TCGAbiolinks`, `scRNAseq`, `biomaRt`, `GseaVis`, `maftools`, `copykat`, `future`, `tidyverse`, `dplyr`, `tidyr`, `tibble`, `readr`, `purrr`, `stringr`, `reshape2`, `openxlsx`, `mgcv`, `ggplot2`, `ggpubr`, `ggrepel`, `ComplexHeatmap`, `patchwork`, `gridExtra` and `viridis`. Exact versions for these packages have not been documented in this README. Add them if original environment records become available.

### Original ssGSEA call

The author provided the following original call, which uses the legacy GSVA interface:

```r
gsva_matrix <- gsva(
  as.matrix(tpmkras1),
  sigList,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)
```

For the ssGSEA method, the effective default weighting exponent is `tau = 0.25` and score normalization is enabled (`ssgsea.norm = TRUE`). The supplied `kcdf` and `abs.ranking` arguments do not affect the ssGSEA calculation; these settings apply to the GSVA algorithm when `method = "gsva"`.

GSVA 1.48 supports this legacy interface. GSVA 1.50 introduced parameter objects while retaining the legacy interface with a deprecation warning; from 1.52, the legacy interface became defunct. An unmodified legacy script therefore requires a compatible GSVA version. Installing the latest GSVA package is not sufficient to reproduce the original call without adapting it.

### Immune-hot/immune-cold clustering

The original clustering call supplied by the author was:

```r
p <- pheatmap(
  gsva_matrix,
  scale = "row",
  show_colnames = TRUE,
  cluster_rows = FALSE,
  fontsize_col = 5
)
cluster <- cutree(p$tree_col, k = 2)
```

Sample clustering was performed on row-scaled ssGSEA scores, and the column dendrogram was cut into two groups. Numeric `cutree()` labels do not themselves identify immune-hot or immune-cold groups; the biological labels should follow the immune-score patterns and the supplied sample metadata. The compatibility of this call with pheatmap 1.0.12 does not establish that 1.0.12 was the original installed version.

Use the original analysis environment where available. Serialized R objects may also depend on package versions. No complete dependency lockfile or independently tested package-version combination is established by the version list above.

## Citation and contact

Please cite the associated manuscript when publicly available, together with the original datasets/methods when using these analyses.
