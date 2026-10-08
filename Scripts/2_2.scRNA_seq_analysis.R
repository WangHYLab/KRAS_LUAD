# ============================================================================
# PART 1 - PACKAGES AND INPUT DATA
# ============================================================================
library(future)
library(Seurat)
library(org.Hs.eg.db)
library(biomaRt)
library(copykat)

# Load the GSE136246 samples
GSM4043245_t1 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043245_NSC015.t1.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043245_t1 <- t(GSM4043245_t1)
GSM4043245_t1obj <- CreateSeuratObject(counts = GSM4043245_t1, project = "GSM4043245_t1", min.cells = 3, min.features = 200)
saveRDS(GSM4043245_t1obj,paste0('/xx /KRAS_LUAD_sc/result/rds_result/filter_GSM4043245_t1.rds'))

GSM4043246_t2 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043246_NSC015.t2.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043246_t2 <- t(GSM4043246_t2)
GSM4043246_t2obj <- CreateSeuratObject(counts = GSM4043246_t2, project = "GSM4043246_t2", min.cells = 3, min.features = 200)
saveRDS(GSM4043246_t2obj,paste0('/xx /KRAS_LUAD_sc/result/rds_result/filter_GSM4043246_t2.rds'))

GSM4043247_t3 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043247_NSC015.t3.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043247_t3 <- t(GSM4043247_t3)
GSM4043247_t3obj <- CreateSeuratObject(counts = GSM4043247_t3, project = "GSM4043247_t3", min.cells = 3, min.features = 200)
saveRDS(GSM4043247_t3obj,paste0('/xx /KRAS_LUAD_sc/result/rds_result/filter_GSM4043247_t3.rds'))

GSM4043248_t1 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043248_NSC016.t1.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043248_t1 <- t(GSM4043248_t1)
GSM4043248_t1obj <- CreateSeuratObject(counts = GSM4043248_t1, project = "GSM4043248_t1", min.cells = 3, min.features = 200)
saveRDS(GSM4043248_t1obj,paste0('/xx /KRAS_LUAD_sc/result/rds_result/filter_GSM4043248_t1.rds'))

GSM4043249_t2 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043249_NSC016.t2.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043249_t2 <- t(GSM4043249_t2)
GSM4043249_t2obj <- CreateSeuratObject(counts = GSM4043249_t2, project = "GSM4043249_t2", min.cells = 3, min.features = 200)
saveRDS(GSM4043249_t2obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043249_t2.rds'))

GSM4043250_t3 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043250_NSC016.t3.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043250_t3 <- t(GSM4043250_t3)
GSM4043250_t3obj <- CreateSeuratObject(counts = GSM4043250_t3, project = "GSM4043250_t3", min.cells = 3, min.features = 200)
saveRDS(GSM4043250_t3obj,paste0('/xx /KRAS_LUAD_sc/result/rds_result/filter_GSM4043250_t3.rds'))

GSM4043251_t1 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043251_NSC019.t1.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043251_t1 <- t(GSM4043251_t1)
GSM4043251_t1obj <- CreateSeuratObject(counts = GSM4043251_t1, project = "GSM4043251_t1", min.cells = 3, min.features = 200)
saveRDS(GSM4043251_t1obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043251_t1.rds'))

GSM4043252_t2 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043252_NSC019.t2.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043252_t2 <- t(GSM4043252_t2)
GSM4043252_t2obj <- CreateSeuratObject(counts = GSM4043252_t2, project = "GSM4043252_t2", min.cells = 3, min.features = 200)
saveRDS(GSM4043252_t2obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043252_t2.rds'))

GSM4043253_t1 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043253_NSC020.t1.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043253_t1 <- t(GSM4043253_t1)
GSM4043253_t1obj <- CreateSeuratObject(counts = GSM4043253_t1, project = "GSM4043253_t1", min.cells = 3, min.features = 200)
saveRDS(GSM4043253_t1obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043253_t1.rds'))

GSM4043254_t2 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043254_NSC020.t2.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043254_t2 <- t(GSM4043254_t2)
GSM4043254_t2obj <- CreateSeuratObject(counts = GSM4043254_t2, project = "GSM4043254_t2", min.cells = 3, min.features = 200)
saveRDS(GSM4043254_t2obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043254_t2.rds'))

GSM4043255_t1 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043255_NSC021.t1.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043255_t1 <- t(GSM4043255_t1)
GSM4043255_t1obj <- CreateSeuratObject(counts = GSM4043255_t1, project = "GSM4043255_t1", min.cells = 3, min.features = 200)
saveRDS(GSM4043255_t1obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043255_t1.rds'))

GSM4043256_t2 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043256_NSC021.t2.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043256_t2 <- t(GSM4043256_t2)
GSM4043256_t2obj <- CreateSeuratObject(counts = GSM4043256_t2, project = "GSM4043256_t2", min.cells = 3, min.features = 200)
saveRDS(GSM4043256_t2obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043256_t2.rds'))

GSM4043257 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043257_NSC035.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043257 <- t(GSM4043257)
GSM4043257obj <- CreateSeuratObject(counts = GSM4043257, project = "GSM4043257", min.cells = 3, min.features = 200)
saveRDS(GSM4043257obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043257.rds'))

GSM4043258 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043258_NSC036.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043258 <- t(GSM4043258)
GSM4043258obj <- CreateSeuratObject(counts = GSM4043258, project = "GSM4043258", min.cells = 3, min.features = 200)
saveRDS(GSM4043258obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043258.rds'))

GSM4043260 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE136246_KRASmut_Gene_vs_Cell raw counts matrix (single cell matrix)/GSM4043260_NSC040.counts.tsv", header = T, row.names = 1, check.names=FALSE)
GSM4043260 <- t(GSM4043260)
GSM4043260obj <- CreateSeuratObject(counts = GSM4043260, project = "GSM4043260", min.cells = 3, min.features = 200)
saveRDS(GSM4043260obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4043260.rds'))

# Load the GSE123902 samples
GSM3516669 <- read.csv("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE123902_KRASmut_ csv of cells x genes unnormalized/GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv", header = T, row.names = 1, check.names=FALSE)
GSM3516669 <- t(GSM3516669)
GSM3516669obj <- CreateSeuratObject(counts = GSM3516669, project = "GSM3516669", min.cells = 3, min.features = 200)
saveRDS(GSM3516669obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM3516669.rds'))

GSM3516672 <- read.csv("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE123902_KRASmut_ csv of cells x genes unnormalized/GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv", header = T, row.names = 1, check.names=FALSE)
GSM3516672 <- t(GSM3516672)
GSM3516672obj <- CreateSeuratObject(counts = GSM3516672, project = "GSM3516672", min.cells = 3, min.features = 200)
saveRDS(GSM3516672obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM3516672.rds'))

GSM3516674 <- read.csv("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE123902_KRASmut_ csv of cells x genes unnormalized/GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv", header = T, row.names = 1, check.names=FALSE)
GSM3516674 <- t(GSM3516674)
GSM3516674obj <- CreateSeuratObject(counts = GSM3516674, project = "GSM3516674", min.cells = 3, min.features = 200)
saveRDS(GSM3516674obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM3516674.rds'))

# Load the normalized GSE146100 sample
GSE146100 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE146100_KRAS(GSM436353 is KRASmut)_TXT file contain counts for each sample/GSE146100_NormData.txt", header = T, row.names = 1, check.names=FALSE)
GSE146100_W2 <- GSE146100[, grep("^W2", names(GSE146100))]
GSM436353obj <- CreateSeuratObject(counts = GSE146100_W2, project = "GSE146100_W2", min.cells = 3, min.features = 200)
saveRDS(GSM436353obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM436353.rds'))

# Load the GSE148466 sample
GSM4472055 <- read.table("/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE148466_KRASmut_Tab-delimited text file. Gene expression (as counts) of each gene for each extracted cell/GSM4472055_P8126999.dge.txt", header = T, row.names = 1, check.names=FALSE)
GSM4472055obj <- CreateSeuratObject(counts = GSM4472055, project = "GSM4472055", min.cells = 3, min.features = 200)
saveRDS(GSM4472055obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4472055.rds'))

# Load the GSE149655 samples
GSM4506699 <- Read10X(data.dir='/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE149655_KRASmut_batched adjusted data of all samples/GSM4506699', gene.column = 1)
GSM4506699 <- Read10X(data.dir='/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE149655_KRASmut_batched adjusted data of all samples/GSM4506699', gene.column = 1)
GSM4506699obj <- CreateSeuratObject(counts = GSM4506699, project ="GSM4506699", min.cells = 3, min.features = 200)
saveRDS(GSM4506699obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4506699.rds'))

GSM4506701 <- Read10X(data.dir='/xx/KRAS_LUAD_sc/KRASmut_scrna_data/GSE149655_KRASmut_batched adjusted data of all samples/GSM4506701', gene.column = 1)
GSM4506701obj <- CreateSeuratObject(counts = GSM4506701, project ="GSM4506701", min.cells = 3, min.features = 200)
saveRDS(GSM4506701obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_GSM4506701.rds'))

# ============================================================================
# PART 2 - QUALITY CONTROL AND PREPROCESSING
# ============================================================================
# Define the sample sets
samples <- c("GSM4043245_t1","GSM4043246_t2","GSM4043247_t3","GSM4043248_t1","GSM4043249_t2","GSM4043250_t3","GSM4043251_t1","GSM4043252_t2","GSM4043253_t1","GSM4043254_t2",
             "GSM4043255_t1","GSM4043256_t2","GSM4043257","GSM4043258","GSM4043260","GSM3516669","GSM3516672", "GSM3516674","GSM436353","GSM4472055", "GSM4506699","GSM4506701")
samples_without_norm <- c("GSM4043245_t1","GSM4043246_t2","GSM4043247_t3","GSM4043248_t1","GSM4043249_t2","GSM4043250_t3","GSM4043251_t1","GSM4043252_t2","GSM4043253_t1","GSM4043254_t2",
             "GSM4043255_t1","GSM4043256_t2","GSM4043257","GSM4043258","GSM4043260","GSM3516669","GSM3516672", "GSM3516674","GSM4472055", "GSM4506699","GSM4506701")

# Apply quality control, normalization, and variable-feature selection
for (sample in samples_without_norm){
  filter_seuratobj <- readRDS(paste0('/xx/KRAS_LUAD_sc/result/rds_result/filter_',sample,'.rds'))
  filter_seuratobj[["percent.mt"]] <- PercentageFeatureSet(filter_seuratobj, pattern = "^MT-")
  filter_seuratobj[["percent.ribo"]] <- PercentageFeatureSet(filter_seuratobj, pattern = "^RP[SL]")
  seuratobj <- subset(filter_seuratobj, subset = nFeature_RNA > 200 & nCount_RNA > 500 &
                        nFeature_RNA < 5000 & nCount_RNA < 20000 &
                        percent.mt < 10 & percent.ribo < 50)
  seuratobj <- NormalizeData(seuratobj)
  seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
  saveRDS(seuratobj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_',sample,'.rds'))
}

# Process the pre-normalized GSE146100 sample
GSM436353obj[["percent.mt"]] <- PercentageFeatureSet(GSM436353obj, pattern = "^MT-")
GSM436353obj[["percent.ribo"]] <- PercentageFeatureSet(GSM436353obj, pattern = "^RP[SL]")
GSM436353obj <- FindVariableFeatures(GSM436353obj, selection.method = "vst", nfeatures = 2000)
saveRDS(GSM436353obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM436353.rds'))

# ============================================================================
# PART 3 - GENE IDENTIFIER CONVERSION
# ============================================================================
# Convert Ensembl identifiers with biomaRt
RenameGenesSeurat <- function(obj ,
                              newnames ) {

  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = mart)
GSM4506699obj <- readRDS('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506699.rds')
gene.ids <- rownames(GSM4506699obj)
gene.names <- genes$external_gene_name[match(gene.ids, genes$ensembl_gene_id)]
replace_indices <- gene.names == '' | is.na(gene.names)
gene.names[replace_indices] <- gene.ids[replace_indices]
GSM4506699obj <- RenameGenesSeurat(obj = GSM4506699obj,
                                   newnames = gene.names)
saveRDS(GSM4506699obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506699.rds'))

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = mart)
GSM4506701obj <- readRDS('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506701.rds')
gene.ids <- rownames(GSM4506701obj)
gene.names <- genes$external_gene_name[match(gene.ids, genes$ensembl_gene_id)]
replace_indices <- gene.names == '' | is.na(gene.names)
gene.names[replace_indices] <- gene.ids[replace_indices]
GSM4506701obj <- RenameGenesSeurat(obj = GSM4506701obj,
                                   newnames = gene.names)
saveRDS(GSM4506701obj,paste0('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506701.rds'))

# Complete the conversion with the org.Hs.eg.db mapping
convert_ensembl_to_symbol <- function(seurat_obj) {
  gene.ids <- rownames(seurat_obj)
  gene.map <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = gene.ids,
                                    columns = "SYMBOL",
                                    keytype = "ENSEMBL")
  gene.names <- gene.map$SYMBOL[match(gene.ids, gene.map$ENSEMBL)]
  gene.names[is.na(gene.names) | gene.names == ""] <- gene.ids[is.na(gene.names) | gene.names == ""]

  RenameGenesSeurat <- function(obj , newnames ) {
    RNA <- obj@assays$RNA
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]         <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]]             <- newnames
      if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
    } else {
      stop("Mismatch in gene number and name vector")
    }
    obj@assays$RNA <- RNA
    return(obj)
  }

  seurat_obj <- RenameGenesSeurat(seurat_obj, gene.names)
  return(seurat_obj)
}

GSM4506699obj <- readRDS("/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506699.rds")
GSM4506699obj <- convert_ensembl_to_symbol(GSM4506699obj)
saveRDS(GSM4506699obj, "/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506699.rds")

GSM4506701obj <- readRDS("/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506701.rds")
GSM4506701obj <- convert_ensembl_to_symbol(GSM4506701obj)
saveRDS(GSM4506701obj, "/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_GSM4506701.rds")

# ============================================================================
# PART 4 - COPYKAT INFERENCE AND METADATA
# ============================================================================
# Run CopyKAT for each sample and attach the predictions to the Seurat objects
setwd("/xx/KRAS_LUAD_sc/result/copykat")
for (sample in samples){
  preprocessed_seuratobj <- readRDS(paste0('/xx/KRAS_LUAD_sc/result/rds_result/preprocessed_strict_',sample,'.rds'))
  exp.rawdata <- GetAssayData(preprocessed_seuratobj, slot = "counts")
  copykat.result <- copykat(exp.rawdata,
                            id.type = "S",
                            ngene.chr = 3,
                            win.size = 25,
                            KS.cut = 0.1,
                            sam.name = sample,
                            distance = "euclidean",
                            n.cores = 4,
                            norm.cell.names = "",
                            output.seg = FALSE)
  mallignant <- read.delim(paste0('/xx/KRAS_LUAD_sc/result/copykat/',sample,'_copykat_prediction.txt'),row.names = 1)
  preprocessed_cnv_seuratobj <- AddMetaData(preprocessed_seuratobj, metadata = mallignant)
  saveRDS(preprocessed_cnv_seuratobj,paste0('/xx/KRAS_LUAD_sc/result0408/rds_result/rds_result_cnv/preprocessed_cnv_',sample,'.rds'))
}

# ============================================================================
# PART 5 - DATA INTEGRATION, CLUSTERING, AND SAMPLE METADATA
# ============================================================================
# Integrate the CopyKAT-annotated samples
preprocessed_seuratobj_list <- list()
for (sample in samples){
  preprocessed_seuratobj <- readRDS(paste0('/xx/KRAS_LUAD_sc/result0408/rds_result/rds_result_cnv/preprocessed_cnv_',sample,'.rds'))
  preprocessed_seuratobj_list[[sample]] <- preprocessed_seuratobj
}
features <- SelectIntegrationFeatures(object.list = preprocessed_seuratobj_list)
plan("multiprocess", workers = 15)
options(future.globals.maxSize = 10000 * 1024^2)
dims=1:20
anchors <- FindIntegrationAnchors(object.list = preprocessed_seuratobj_list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors,dims = dims)

# Scale the integrated assay and calculate PCA
DefaultAssay(combined) <- "integrated"
dims=1:20
seurat.combined <- ScaleData(combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined,  features = VariableFeatures(object = seurat.combined))

# Build the neighbor graph, identify clusters, and calculate UMAP
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:20)
seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)
seurat.combined <- RunUMAP(seurat.combined, dims = 1:20)

# Preserve the RNA assay as the default assay in the final object
DefaultAssay(seurat.combined) <- "RNA"

# Standardize sample identifiers and their display order
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043245_t1"), "sample"] = "GSM4043245_t1"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043246_t2"), "sample"] = "GSM4043246_t2"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043247_t3"), "sample"] = "GSM4043247_t3"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043248_t1"), "sample"] = "GSM4043248_t1"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043249_t2"), "sample"] = "GSM4043249_t2"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043250_t3"), "sample"] = "GSM4043250_t3"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043251_t1"), "sample"] = "GSM4043251_t1"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043252_t2"), "sample"] = "GSM4043252_t2"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043253_t1"), "sample"] = "GSM4043253_t1"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043254_t2"), "sample"] = "GSM4043254_t2"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043255_t1"), "sample"] = "GSM4043255_t1"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043256_t2"), "sample"] = "GSM4043256_t2"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043257"), "sample"] = "GSM4043257"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043258"), "sample"] = "GSM4043258"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4043260"), "sample"] = "GSM4043260"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM3516669"), "sample"] = "GSM3516669"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM3516672"), "sample"] = "GSM3516672"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM3516674"), "sample"] = "GSM3516674"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="W2"), "sample"] = "GSM436353"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4472055"), "sample"] = "GSM4472055"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4506699"), "sample"] = "GSM4506699"
seurat.combined@meta.data[which(seurat.combined@meta.data$orig.ident=="GSM4506701"), "sample"] = "GSM4506701"

seurat.combined@meta.data$sample <- factor(seurat.combined@meta.data$sample,
                                                     levels=c("GSM4043245_t1","GSM4043246_t2","GSM4043247_t3",
                                                              "GSM4043248_t1","GSM4043249_t2","GSM4043250_t3",
                                                              "GSM4043251_t1","GSM4043252_t2",
                                                              "GSM4043253_t1","GSM4043254_t2",
                                                              "GSM4043255_t1","GSM4043256_t2",
                                                              "GSM4043257","GSM4043258","GSM4043260",
                                                              "GSM3516669","GSM3516672", "GSM3516674",
                                                              "GSM436353","GSM4472055", "GSM4506699","GSM4506701") )

# ============================================================================
# PART 6 - FINAL CELL-TYPE ANNOTATION AND OUTPUT
# ============================================================================
# Apply the established broad cell-type identities
seurat.combined2 <- seurat.combined
personal.cluster.ids <- c("immun","immun","immun","immun","epithelial","immun","immun","immun","immun","immun","immun","immun","stromal","stromal","immun","immun","epithelial")
names(personal.cluster.ids) <- levels(seurat.combined2)
seurat.combined2 <- RenameIdents(seurat.combined2, personal.cluster.ids)

# Create the final object and annotate malignant, immune, and stromal cells
seurat.combined9 <- seurat.combined2
seurat.combined9@meta.data$copykat_seurat_cluster <- NA
seurat.combined9@meta.data[which(seurat.combined9@meta.data$seurat_clusters %in% c('4', '16')),"copykat_seurat_cluster"] <- "malignant"
seurat.combined9@meta.data[which(seurat.combined9@meta.data$seurat_clusters %in%
                                   c('0', '1', '2', '3', '5', '6', '7', '8', '9', '10', '11', '14', '15')),
                           "copykat_seurat_cluster"] <- "immun"
seurat.combined9@meta.data[which(seurat.combined9@meta.data$seurat_clusters %in%
                                   c('12', '13')),
                           "copykat_seurat_cluster"] <- "stromal"
saveRDS(seurat.combined9, file = "/xx/KRAS_LUAD_sc/result/rds_result2/seurat_combined9_updated.rds")

# ============================================================================
# PART 7 - SCALED MALIGNANT-CELL EXPRESSION AND CORRELATION ANALYSIS
# ============================================================================
# Create a temporary analysis sample field without changing the saved final RDS
analysis_sample <- as.character(seurat.combined9@meta.data$sample)
gse146100_cells <- seurat.combined9@meta.data$orig.ident %in% c("W2", "GSE146100_W2")
analysis_sample[gse146100_cells] <- "GSM436353"

seurat.analysis <- seurat.combined9
seurat.analysis@meta.data$analysis_sample <- factor(
  analysis_sample,
  levels = samples
)
scRNAlist <- SplitObject(seurat.analysis, split.by = "analysis_sample")

load("/xx/KRAS_LUAD_sc/data/exp_count502.RData")
gene.all <- unique(exp_count502$gene_name)

immune_percentages <- setNames(rep(NA_real_, length(samples)), samples)
malignant_expression_means <- matrix(
  NA_real_,
  nrow = length(gene.all),
  ncol = length(samples),
  dimnames = list(gene.all, samples)
)

# Scale each sample separately and summarize expression in malignant cells
for (sample in samples) {
  sample_obj <- scRNAlist[[sample]]
  if (is.null(sample_obj)) {
    next
  }

  sample_metadata <- sample_obj@meta.data
  immune_percentages[sample] <- sum(
    sample_metadata$copykat_seurat_cluster == "immun",
    na.rm = TRUE
  ) / nrow(sample_metadata)

  original_obj <- readRDS(
    paste0(
      "/xx/KRAS_LUAD_sc/result/rds_result/",
      "preprocessed_strict_",
      sample,
      ".rds"
    )
  )
  features_to_scale <- Reduce(
    intersect,
    list(
      rownames(original_obj),
      rownames(sample_obj),
      gene.all
    )
  )

  sample_obj <- ScaleData(
    sample_obj,
    assay = "RNA",
    features = features_to_scale,
    verbose = FALSE
  )

  malignant_cells <- rownames(
    sample_metadata[
      sample_metadata$copykat_seurat_cluster == "malignant",
      ,
      drop = FALSE
    ]
  )
  malignant_cells <- intersect(
    malignant_cells,
    colnames(sample_obj@assays$RNA@scale.data)
  )

  if (length(malignant_cells) > 0) {
    malignant_expression_means[features_to_scale, sample] <- rowMeans(
      sample_obj@assays$RNA@scale.data[
        features_to_scale,
        malignant_cells,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
}

# Calculate the across-sample Pearson correlation for every gene
correlation_r <- rep(NA_real_, length(gene.all))
correlation_p <- rep(NA_real_, length(gene.all))
malignant_expression_text <- character(length(gene.all))
immune_percentage_text <- paste(immune_percentages, collapse = ",")

for (i in seq_along(gene.all)) {
  gene_values <- malignant_expression_means[i, ]
  malignant_expression_text[i] <- paste(gene_values, collapse = ",")
  valid <- complete.cases(immune_percentages, gene_values)

  if (
    sum(valid) > 2 &&
    length(unique(immune_percentages[valid])) > 1 &&
    length(unique(gene_values[valid])) > 1
  ) {
    correlation_result <- cor.test(
      immune_percentages[valid],
      gene_values[valid],
      method = "pearson"
    )
    correlation_r[i] <- unname(correlation_result$estimate)
    correlation_p[i] <- correlation_result$p.value
  }
}

unique_df <- data.frame(
  r = signif(correlation_r, digits = 4),
  p = signif(correlation_p, digits = 4),
  immune_percentage = rep(immune_percentage_text, length(gene.all)),
  malignant_cell_exp_mean = malignant_expression_text,
  gene_name = gene.all,
  stringsAsFactors = FALSE
)
unique_df$p.adj <- p.adjust(unique_df$p, method = "BH")
unique_df <- unique_df[order(unique_df$p, na.last = TRUE), , drop = FALSE]
rownames(unique_df) <- unique_df$gene_name

write.csv(
  unique_df,
  file = "/xx/KRAS_LUAD_sc/result/statistic/sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv",
  row.names = TRUE
)

# ============================================================================
# PART 8 - FINAL UMAP PLOTS AND COORDINATE EXPORT
# ============================================================================
# Add a plotting-only sample field; the saved final RDS remains unchanged
seurat.combined9@meta.data$sample2 <- factor(
  analysis_sample,
  levels = samples
)

p_umap_sample2 <- DimPlot(
  seurat.combined9,
  reduction = "umap",
  group.by = "sample2"
)
p_umap_celltype <- DimPlot(
  seurat.combined9,
  reduction = "umap",
  group.by = "copykat_seurat_cluster"
)

p_umap_sample2
p_umap_celltype

umap_df <- as.data.frame(
  Embeddings(seurat.combined9, reduction = "umap")
)
umap_df$cell <- rownames(umap_df)
umap_df$sample2 <- seurat.combined9@meta.data[
  rownames(umap_df),
  "sample2"
]
umap_df$celltype <- seurat.combined9@meta.data[
  rownames(umap_df),
  "copykat_seurat_cluster"
]

write.csv(
  umap_df,
  file = "/xx/KRAS_LUAD_sc/result/statistic/umap_df_samples.csv",
  row.names = FALSE
)
