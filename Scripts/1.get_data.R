# =============================================================================
# 1. Set up the working directory and load packages
# =============================================================================
setwd('/path/to/get_data')
#BiocManager::install("TCGAbiolinks",force = TRUE)
#BiocManager::install("GenomeInfoDb",force = TRUE)
#BiocManager::install("scRNAseq",force = TRUE)
#BiocManager::install("Biostrings",force = TRUE)
library(TCGAbiolinks)
library(scRNAseq)

# =============================================================================
# 2. Download TCGA-LUAD gene expression data
# =============================================================================
LUAD<-GDCquery(project="TCGA-LUAD",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "STAR - Counts")
GDCdownload(LUAD,method="api")
exp<-GDCprepare(query=LUAD)
save(exp,file='exp.RData')
load(file='exp.RData')
data<-as.data.frame(rowRanges(exp))

# =============================================================================
# 3. Extract expression matrices and select expression samples
# =============================================================================
# Extract counts with gene annotations.
counts<-as.data.frame(assay(exp))
exp_count<-cbind(genetype=data$gene_type,gene_name=data$gene_name,counts)
save(exp_count,file='exp_count.RData')
load(file='exp_count.RData')
# Find duplicate 16-character sample IDs.
exps16<-substr(colnames(exp_count[3:602]),1,16)
dupexps16<-unique(exps16[duplicated(exps16)])
select_samples_exp<-c()
undups_exp<-c()
for(i in colnames(exp_count[3:602])){
  if(substr(i,1,16) %in% dupexps16){
    select_samples_exp<-c(select_samples_exp,i)
  } else{
    undups_exp<-c(undups_exp,i)
  }
}
sse<-sort(select_samples_exp)
# Keep the selected even-indexed duplicate barcodes.
sample1<-sse[c(2,4,6,8,10,12,14,16,18,20,22)]
samples_exp<-c(undups_exp,sample1)
exp_count_dedup<-exp_count[,c("genetype","gene_name",samples_exp)]
save(exp_count_dedup,file='exp_count_dedup.RData')
load(file='exp_count_dedup.RData')
# Extract unstranded TPM using the same sample selection.
tpm<-as.data.frame(assay(exp,i="tpm_unstrand"))
exp_tpm<-cbind(genetype=data$gene_type,gene_name=data$gene_name,tpm)
save(exp_tpm,file='exp_tpm.RData')
load(file='exp_tpm.RData')
exp_tpm_dedup<-exp_tpm[,c("genetype","gene_name",samples_exp)]
save(exp_tpm_dedup,file='exp_tpm_dedup.RData')
load(file='exp_tpm_dedup.RData')


# =============================================================================
# 4. Process mutation data and export mutation-matched expression cohorts
# =============================================================================

# =============================================================================
# 4.1. Download TCGA-LUAD somatic mutation data
# =============================================================================
LUAD<-GDCquery(project="TCGA-LUAD",data.category = "Simple Nucleotide Variation",data.type = "Masked Somatic Mutation",access = "open")
GDCdownload(LUAD,method="api")
mute<-GDCprepare(query=LUAD)
save(mute,file='mute.RData')

# =============================================================================
# 4.2. Construct a gene-by-sample mutation classification matrix
# =============================================================================
load(file='mute.RData')
genes <- unique(mute$Hugo_Symbol)
samples <- unique(mute$Tumor_Sample_Barcode)
# Zero indicates no recorded mutation.
mutation_matrix <- matrix(0, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))
# Store all mutation classes; the last record wins for each gene-sample pair.
for (i in 1:nrow(mute)) {
  gene <- mute$Hugo_Symbol[i]
  sample <- mute$Tumor_Sample_Barcode[i]
  mutation_matrix[gene, sample] <- mute$Variant_Classification[i]
}
save(mutation_matrix,file='mutation_matrix.RData')
load(file='mutation_matrix.RData')
mutation_matrix<-as.data.frame(mutation_matrix)
# =============================================================================
# 4.3. Select mutation records using 16-character sample identifiers
# =============================================================================
# Find duplicate 16-character sample IDs.
mutes16<-substr(colnames(mutation_matrix),1,16)
dupmutes16<-unique(mutes16[duplicated(mutes16)])
select_samples<-c()
undups<-c()
for(i in colnames(mutation_matrix)){
  if(substr(i,1,16) %in% dupmutes16){
    select_samples<-c(select_samples,i)
  } else{
    undups<-c(undups,i)
  }
}
# Export KRAS/EGFR classes for duplicate-sample review.
muted<-t(mutation_matrix[,select_samples])
muted_k<-as.data.frame(muted[,c('KRAS','EGFR')])
mutedata<-muted_k[order(rownames(muted_k)),]
write.csv(mutedata,file='mutedata.csv')

# Manually retain one barcode per duplicate sample in mutedata_modified.csv (X).
mutedm<-read.csv('mutedata_modified.csv')
mutes21<-mutedm$X
mute_samples_F<-c(undups,mutes21)

mutation_matrix_dedup<-mutation_matrix[,mute_samples_F]
save(mutation_matrix_dedup,file='mutation_matrix_dedup.RData')

# Select samples with a recorded KRAS mutation.
mkras<-as.data.frame(t(mutation_matrix_dedup["KRAS",]))
mkras1<-mkras[which(mkras$KRAS!=0),,drop = FALSE]

# =============================================================================
# 4.4. Match mutation and expression samples and assign KRAS status
# =============================================================================
# Match mutation and expression samples by 16-character IDs.
mute_samples<-substr(colnames(mutation_matrix_dedup),1,16)
exp_samples<-substr(colnames(exp_count_dedup)[3:591],1,16)
m_e_s<-intersect(mute_samples,exp_samples)
k<-intersect(substr(rownames(mkras1),1,16),m_e_s)
wt<-setdiff(m_e_s,k)
# Label KRAS-mutant samples 1 and others 0 (fixed sizes: 142/377).
samples_KRAS<-data.frame(sample=c(k,wt),KRAS=c(rep('1',142),rep('0',377)))

# Keep primary tumors with vial A (01A).
samples_KRAS1<-samples_KRAS[which(substr(samples_KRAS$sample,16,16)=='A'),]
sk<-samples_KRAS1[which(substr(samples_KRAS1$sample,14,16)=='01A'),]
write.csv(sk,file='sk.csv',row.names = F)

# =============================================================================
# 4.5. Export expression datasets for the matched primary tumor cohort
# =============================================================================
# Shorten barcodes to 16 characters and align samples to sk$sample.
colnames(exp_count_dedup)<-c("genetype","gene_name",substr(colnames(exp_count_dedup[3:591]),1,16))
exp_count502<-exp_count_dedup[,c("genetype","gene_name",sk$sample)]
save(exp_count502,file='exp_count502.RData')
colnames(exp_tpm_dedup)<-c("genetype","gene_name",substr(colnames(exp_tpm_dedup[3:591]),1,16))
exp_tpm502<-exp_tpm_dedup[,c("genetype","gene_name",sk$sample)]
save(exp_tpm502,file='exp_tpm502.RData')

# =============================================================================
# 4.6. Export expression datasets for the KRAS-mutant subset
# =============================================================================
k<-sk[which(sk$KRAS==1),]
kras_count138<-exp_count_dedup[,c("genetype","gene_name",k$sample)]
kras_tpm138<-exp_tpm_dedup[,c("genetype","gene_name",k$sample)]
save(kras_count138,file='kras_count138.RData')
save(kras_tpm138,file='kras_tpm138.RData')

# =============================================================================
# 4.7. Export mutation classifications for the matched primary tumor cohort
# =============================================================================
# Align mutation samples to the expression cohort.
colnames(mutation_matrix_dedup)<-c(substr(colnames(mutation_matrix_dedup),1,16))
mutation502<-mutation_matrix_dedup[,sk$sample]
save(mutation502,file='mutation502.RData')

