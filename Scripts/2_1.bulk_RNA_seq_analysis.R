setwd('path/to/data')
WORK_DIR<-"path/to/KRAS_LUAD/"
# TCGA-LUAD bulk RNA-seq analysis: KRAS mutation, immune enrichment,
# immune subtyping, subtype comparisons, and TIMER2 visualization.
# Figure map:
#   Fig. 2a: CIBERSORT GSEA bubble plot
#   Fig. 2b: ImmuCellAI GSEA bubble plot
#   Fig. 2c: CIBERSORT GSEA curves
#   Fig. 2d: ImmuCellAI GSEA curves
#   Fig. 3a: ssGSEA subtype heatmap
#   Fig. 3b: subtype-associated mutation oncoplot
#   Fig. 3c: immunotherapy-marker density plot
#   Fig. 3d: subtype survival curves
#   Fig. 3e: GO BP enrichment 
#   Supplementary Fig. 1a: TIMER2 log2FC heatmap
# Fixed sample counts and geneSetID values should be checked if inputs change.
library(dplyr)
#BiocManager::install("DESeq2",force = TRUE)
library(DESeq2)
# ==============================================================================
# Part 1. Prepare counts, TPM, and mutation data; remove duplicated gene names
# ==============================================================================
load(file='exp_count502.RData')
load(file='exp_tpm502.RData')
load(file='mutation502.RData')
exp_count<-exp_count502[!duplicated(exp_count502$gene_name),]
exp_count<-exp_count[,-1]
save(exp_count,file='exp_count.RData')
exp_tpm<-exp_tpm502[!duplicated(exp_tpm502$gene_name),]
save(exp_tpm,file='exp_tpm.RData')
exp_tpm_coding<-exp_tpm[which(exp_tpm$genetype=='protein_coding'),]
exp_tpm_c<-exp_tpm_coding[,-1]
save(exp_tpm_c,file='exp_tpm_c.RData')


# ==============================================================================
# Part 2. KRAS-mutant versus wild-type analysis (Fig. 2a-d)
# ==============================================================================
# 2.1. Differential expression: mut versus wt.
# Define KRAS-mut and KRAS-wt samples.
krasm<-mutation502[,which(mutation502['KRAS',]!=0)]
samples<-c(colnames(krasm),setdiff(colnames(mutation502),colnames(krasm)))
type<-factor(c(rep('mut',138),rep('wt',364)))
colData<-data.frame(row.names=samples,type=type)

## Prepare protein-coding count data for DESeq2.
gene_c<-exp_count502[which(exp_count502$genetype=='protein_coding'),]$gene_name
rownames(exp_count)<-exp_count$gene_name
count<-exp_count[,-1]
count1<-subset(count,select=rownames(colData))
count2<-as.matrix(count1[gene_c,])

## DESeq2: positive log2FC indicates higher expression in mut.
library(DESeq2)
dds<-DESeqDataSetFromMatrix(count2, colData, design= ~ type)
dds$type<-relevel(dds$type,ref='wt')
#save(dds,file='dds.RData')
dds2<-DESeq(dds)
# Reload saved results, replacing the fit above.
#save(dds2,file='dds2.RData')
load(file='dds2.RData')
res<-results(dds2)
res<-res[order(res$padj),]
res<-data.frame(gene_name=rownames(res),res)
write.csv(res,file= "all_deseq2_results_mut_vs_wt.csv",row.names = F)

# 2.2. Immune-signature GSEA (Fig. 2a-d).
# Rank genes by the DESeq2 Wald statistic for GSEA.
deseq_res<-read.csv(file= "all_deseq2_results_mut_vs_wt.csv")
sym_genelist<-deseq_res$stat
names(sym_genelist)<-deseq_res$gene_name
sym_genelist<-sort(sym_genelist,decreasing = T)

## CIBERSORT GSEA
library(clusterProfiler)
# Marker file must follow the TERM2GENE format.
immuneSig<-read.csv('CIBERSORT_cell_types_markers.csv')
immune_ges_stat_cibersort<-GSEA(sym_genelist,TERM2GENE=immuneSig,pvalueCutoff = 1)
save(immune_ges_stat_cibersort,file='immune_ges_stat_cibersort.RData')
load(file='immune_ges_stat_cibersort.RData')
result<-immune_ges_stat_cibersort@result
write.csv(result,file='GSEA_mut_wt_result_cibersort_stat_order.csv',row.names = F)
# ==================== Fig. 2a ====================
# CIBERSORT GSEA bubble plot: x = -log10(q value), size = |NES|.
library(ggplot2)
res<-read.csv(file='GSEA_mut_wt_result_cibersort_stat_order.csv')
res1<-res
res1$annotation<-c(rep('NES<0',22))
res1$logq<--log10(res1$qvalue)
p_ciber<-ggplot(res1,aes(x=logq,y=reorder(ID,logq),size=abs(NES),colour=annotation))+
  scale_color_manual(values=c("blue","red"))+
  geom_point(shape=16)+labs(x='-log10(qvalue)',y='Immune cells')+
  theme_set(theme_bw())+theme(axis.title = element_text(size=13),axis.title.x = element_text(size=13),legend.title = element_text(size=9))
# ==================== Fig. 2c ====================
# CIBERSORT GSEA curves for gene sets 1-4.
library(GseaVis)
a<-gseaNb(object = immune_ges_stat_cibersort,geneSetID = 1,pCol = 'black',pHjust = 0)
b<-gseaNb(object = immune_ges_stat_cibersort,geneSetID = 2,pCol = 'black',pHjust = 0)
c<-gseaNb(object = immune_ges_stat_cibersort,geneSetID = 3,pCol = 'black',pHjust = 0)
d<-gseaNb(object = immune_ges_stat_cibersort,geneSetID = 4,pCol = 'black',pHjust = 0)

plot_list(a,b,c,d,nrow=1)


# Prepare ImmuCellAI signatures for GSEA (Fig. 2b, 2d).
# Convert ImmuCellAI signatures to TERM2GENE format.
library(clusterProfiler)
sig<-read.csv(paste0(WORK_DIR,'Tables/ImmuCellAI_24_immune_cell_types_signature_genes.csv'))
rownames(sig)<-sig[,1]
sig<-sig[,-1]
sig<-as.matrix(sig)
sigList<-list()
for(i in 1:nrow(sig)){
  x<-as.vector(sig[i,])
  x<-x[nchar(x)>0]
  sigList[[i]]<-x
}
cell<-rownames(sig)
names(sigList)<-cell

gene_counts <- sapply(sigList, length)
sig_df <- data.frame(
  term = rep(names(sigList), times = gene_counts),
  gene = unlist(sigList, use.names = FALSE)
)
# Rank genes by the Wald statistic and run ImmuCellAI GSEA.
luad<-read.csv(paste0(WORK_DIR,'/Results/all_deseq2_results_mut_vs_wt.csv'))
stat_luad <- luad$stat
names(stat_luad) <- luad$gene_name
stat_luad <- sort(stat_luad, decreasing = TRUE)
immune_ges_luad24_stat<-GSEA(stat_luad,TERM2GENE=sig_df,pvalueCutoff = 1)
save(immune_ges_luad24_stat,file=paste0(WORK_DIR,'Rdata/immune_ges_luad24_stat.Rdata'))

load(file=paste0(WORK_DIR,'Rdata/immune_ges_luad24_stat.Rdata'))
write.csv(immune_ges_luad24_stat,file=paste0(WORK_DIR,'Results/immune_ges_luad24_stat.csv'),row.names = F)


library(ggplot2)
immune_res_luad<-immune_ges_luad24_stat@result
immune_res_luad$annotation<-c(rep('NES<0',20))
immune_res_luad[which(immune_res_luad$NES>0),'annotation']<-'NES>0'
colnames(immune_res_luad)[c(1,4)]<-c('immunecell','ES')
immune_res_luad$logqvalue<--log10(immune_res_luad$qvalue)

immune_res_luad_sorted <- immune_res_luad %>%
  arrange(logqvalue, desc(Description)) %>%
  # Preserve the sorted cell-type order in the plot.
  mutate(immunecell = factor(immunecell, levels = unique(immunecell)))
# ==================== Fig. 2b ====================
# ImmuCellAI GSEA bubble plot: size = |NES|; color = NES direction.
p_immuncellAI<-ggplot(immune_res_luad_sorted,aes(x=logqvalue,y=immunecell,size=abs(NES),colour=annotation))+
  scale_color_manual(values=c("blue","red"))+
  geom_point(shape=16)+
  #labs(x='-log10(qvalue)',y='Immune cells',title = 'KRAS-mut vs KRAS-wt (TCGA-LUAD)')+
  labs(x='-log10(qvalue)',y='Immune cells')+
  theme_set(theme_bw())+theme(axis.title = element_text(size=13),axis.title.x = element_text(size=13),legend.title = element_text(size=9))
library(GseaVis)
library(patchwork)
# ==================== Fig. 2d ====================
# ImmuCellAI GSEA curves for gene sets 1-4.
# Verify that these IDs match monocyte, gamma, B cell, and Tfh.
p_monocyte<-gseaNb(object = immune_ges_luad24_stat,geneSetID = 1,pCol = 'black',pHjust = 0)
p_gamma<-gseaNb(object = immune_ges_luad24_stat,geneSetID = 2,pCol = 'black',pHjust = 0)
p_bcell<-gseaNb(object = immune_ges_luad24_stat,geneSetID = 3,pCol = 'black',pHjust = 0)
p_tfh<-gseaNb(object = immune_ges_luad24_stat,geneSetID = 4,pCol = 'black',pHjust = 0)

plot_list(p_monocyte,p_gamma,p_bcell,p_tfh,nrow=1)

library(gridExtra)
# Combine Fig. 2a and Fig. 2b.
grid.arrange(p_ciber,p_immuncellAI,nrow=1,widths=c(5,4))
# ==============================================================================
# Part 3. Immune subtyping of KRAS-mutant tumors (Fig. 3a)
# ==============================================================================
# Log2-transform TPM values and retain KRAS-mut samples.
tpm<-exp_tpm_c
rownames(tpm)<-tpm$gene_name
tpm1<-tpm[,-1]
tpm1<-log2(tpm1+1)
krasm<-mutation502[,which(mutation502['KRAS',]!=0)]
tpmkras<-subset(tpm1,select=colnames(krasm))
tpmkras1<-tpmkras[which(rowSums(tpmkras)!=0),]
## ImmuCellAI signatures for 24 immune-cell types.
sig<-read.csv('ImmuCellAI_24_immune_cell_types_signature_genes.csv')
rownames(sig)<-sig[,1]
sig<-sig[,-1]
sig<-as.matrix(sig)
sigList<-list()
for(i in 1:nrow(sig)){
  x<-as.vector(sig[i,])
  x<-x[nchar(x)>0]
  sigList[[i]]<-x
}
cell<-rownames(sig)
names(sigList)<-cell

# 3.1. Calculate ssGSEA scores.
# Rows are signatures; columns are KRAS-mut samples.
library(GSVA)
gsva_matrix<-gsva(as.matrix(tpmkras1),sigList,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
save(gsva_matrix,file='gsva_matrix.RData')
write.csv(gsva_matrix,file='gsva_matrix.csv')
# 3.2. Define cold/hot subtypes from sample clustering.
library(pheatmap)
load(file='gsva_matrix.RData')
# ==================== Fig. 3a ====================
# Cluster row-scaled scores; split samples into two subtypes.
p<-pheatmap(gsva_matrix,scale="row",show_colnames = T,cluster_rows = F,fontsize_col=5)
order<-p$tree_col$order
p_samplelist<-colnames(gsva_matrix)[order]
save(p_samplelist,file='p_samplelist.RData')
load(file='p_samplelist.RData')
# Label cluster 1 as cold and cluster 2 as hot.
cluster<-cutree(p$tree_col,k=2)
cluster1<-data.frame(cluster)
cold<-subset(cluster1,cluster=="1")
cold_id<-data.frame(sample=rownames(cold),subtype=rep('cold',66))
hot<-subset(cluster1,cluster=="2")
hot_id<-data.frame(sample=rownames(hot),subtype=rep('hot',72))
kras_sample<-rbind(cold_id,hot_id)
write.csv(kras_sample,file='kras_sample.csv',row.names = F)

# Redraw Fig. 3a with subtype annotation.
load(file='gsva_matrix.RData')
kras_sample<-read.csv(file='kras_sample.csv')
rownames(kras_sample)<-kras_sample$sample
annotation_col<-kras_sample[,-1,drop=F]
p<-pheatmap(gsva_matrix,scale="row",show_colnames = F,cluster_rows = F,fontsize_col=5,annotation_col = annotation_col)
#p<-pheatmap(gsva_matrix,scale="row",show_colnames = F,cluster_rows = T,fontsize_col=5,annotation_col = annotation_col)
order<-p$tree_col$order
p_samplelist<-colnames(gsva_matrix)[order]
save(p_samplelist,file='p_samplelist.RData')

# 3.3. Compare overall survival (Fig. 3d).
# Merge overall-survival data with cold/hot subtype labels.
sur<-read.table('survival_LUAD_survival.txt',sep='\t',header=T)
sur1<-sur[order(sur$sample),]
sur1$sample<-paste(sur1$sample,'A',sep='')
survival502<-sur1[which(sur1$sample %in% colnames(exp_tpm_c)[-1]),]

surv<-survival502
os_data<-cbind(surv$sample,surv$OS,surv$OS.time)
os_data<-as.data.frame(os_data)
colnames(os_data)<-c("sample","OS","OS.time")
osdata<-os_data[which(os_data$sample %in% kras_sample$sample),]
survival_data<-merge(kras_sample,osdata,by='sample')
survival_data$OS<-as.numeric(survival_data$OS)
survival_data$OS.time<-as.numeric(survival_data$OS.time)
#save(survival_data,file='survival_data.RData')
# Reload saved survival data, replacing the table above.
load(file='survival_data.RData')

# ==================== Fig. 3d ====================
# Kaplan-Meier curves with a log-rank P value.
library(survival)
library(ggpubr)
library(survminer)
fit<-survfit(Surv(OS.time,OS)~subtype,data=survival_data)
ggsurvplot(fit, data = survival_data,
           pval = TRUE,
           palette=c("blue", "red"),
           legend.labs=c("cold","hot"),
           legend.title="subtype",
           ylab="Cumulative survival (percentage)",xlab = " Time (Days)",
           conf.int = FALSE
)

# 3.4. Differential expression: cold (66) versus hot (72).
kras_sample<-read.csv('kras_sample.csv')
names<-kras_sample$sample
subtype<-factor(kras_sample$subtype,levels = c('hot','cold'))
colData<-data.frame(row.names=names, subtype=subtype)

# Select protein-coding count data from KRAS-mut samples.
rownames(exp_count)<-exp_count$gene_name
exp_count1<-exp_count[,-1]
kras_count<-subset(exp_count1,select=names)
genes_c<-exp_count502[which(exp_count502$genetype=='protein_coding'),]$gene_name
kras_count<-kras_count[genes_c,]

# DESeq2: positive log2FC indicates higher expression in cold.
library(DESeq2)
dds_ch<-DESeqDataSetFromMatrix(as.matrix(kras_count), colData, design= ~ subtype)
dds_ch$subtype<-relevel(dds_ch$subtype,ref='hot')
save(dds_ch,file='dds_ch.RData')
dds2_ch<-DESeq(dds_ch)
resultsNames(dds2_ch) 
save(dds2_ch,file='dds2_ch.RData')
# Check coefficient names and reference level.
#load(file='dds2_ch.RData')
dds2_ch$subtype
resultsNames(dds2_ch) 
res_ch<-results(dds2_ch)
res_ch<-res_ch[order(res_ch$padj),]
res_ch<-data.frame(gene_name=rownames(res_ch),res_ch)
write.csv(res_ch,file= "all_deseq2_results_cold_vs_hot.csv",row.names = F)


# ==============================================================================
# Part 4. Test mutation-subtype associations and plot selected genes (Fig. 3b)
# ==============================================================================
# Treat Silent as 0 and other nonzero events as mutations.
load('mutation502.RData')
kras_sample<-read.csv('kras_sample.csv')
krasmut<-as.data.frame(t(mutation502[,kras_sample$sample]))
krasmut[krasmut=='Silent']<-'0'
krasmut[krasmut!='0']<-'1'
krasmut1<-data.frame(sample=rownames(krasmut),krasmut)
mutedata<-merge(kras_sample,krasmut1,by='sample')
genes<-colnames(mutedata)[3:17707]
####### Chi-square tests
# Retain chi-square results when the minimum expected count is >=5.
results <- sapply(genes, function(gene) {
  contingency_table <- table(mutedata[[gene]], mutedata$subtype)
  chi_square_result <- chisq.test(contingency_table)
  min_expected_freq <- min(chi_square_result$expected)
  p_value <- chi_square_result$p.value
  return(c(P_Value = p_value, Min_Expected_Frequency = min_expected_freq))
}, simplify = "data.frame")

result_df <- as.data.frame(t(results))
result_df<-result_df[order(result_df$P_Value),]
result_chi<-result_df[which(result_df$Min_Expected_Frequency>=5),]
result_chisq<-data.frame(genename=rownames(result_chi),p_value=result_chi$P_Value,test_method=rep('chisq',nrow(result_chi)))

# Use Fisher's exact test for the remaining genes.
fisher_genes<-setdiff(rownames(result_df),rownames(result_chi))

results_fi <- sapply(fisher_genes, function(gene) {
  contingency_table <- table(mutedata[[gene]], mutedata$subtype)
  fisher_result <- fisher.test(contingency_table)
  p_value <-fisher_result$p.value
  return(p_value)})  
result_fisher<-data.frame(genename=names(results_fi),p_value=results_fi,test_method=rep('fisher',length(results_fi)))
result_mut<-rbind(result_chisq,result_fisher)

############# Count mutated samples in each subtype
cold<-as.data.frame(t(mutedata[which(mutedata$subtype=='cold'),]))
cold1<-cold[-c(1,2),]
cold2<-as.data.frame(lapply(cold1,as.numeric))
rownames(cold2)<-rownames(cold1)
cold3<-rowSums(cold2[rownames(cold2),])
cold_count<-data.frame(genename=rownames(cold2),cold_mut=cold3)

hot<-as.data.frame(t(mutedata[which(mutedata$subtype=='hot'),]))
hot1<-hot[-c(1,2),]
hot2<-as.data.frame(lapply(hot1,as.numeric))
rownames(hot2)<-rownames(hot1)
hot3<-rowSums(hot2[rownames(hot2),])
hot_count<-data.frame(genename=rownames(hot2),hot_mut=hot3)
sample_count<-merge(cold_count,hot_count,by='genename')

res_mutation_cor<-merge(result_mut,sample_count,by='genename')
res_mutation_cor<-res_mutation_cor[order(res_mutation_cor$p_value),]
write.csv(res_mutation_cor,file='res_mutation_corralation.csv',row.names = F)

# Prepare selected mutations for Fig. 3b.
# Exclude Silent events; label other recorded events as mutation.
load(file='p_samplelist.RData')
load(file="mutation502.RData")
library(maftools)
library(ComplexHeatmap)
library(reshape2)
krasmut<-mutation502[,kras_sample$sample]
genes<-c("TPTE","SLC39A12","KEAP1","GRIK2","MMP16","STK11","AFF2","KCNU1","SPEG")
mutation<-krasmut[which(rownames(krasmut) %in% genes),]
mutation[mutation=='Silent']<-'0'
mutation[mutation!='0']<-'mutation'
m<-as.data.frame(t(mutation))
m$sample<-rownames(m)
m1<-merge(m,kras_sample,by='sample')
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = "white", col = NA))
  },
  mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "black", col = NA))
  })
col <- c(Missense_Mutation = "black", Nonsense_Mutation = "blue",Splice_Region='red',Splice_Site='green')
# Order samples by mutation status, then by subtype.
m2<-m1[order(m1$STK11,m1$MMP16,m1$KEAP1,m1$AFF2,m1$TPTE,m1$KCNU1,m1$SPEG,m1$SLC39A12,m$GRIK2,decreasing=T),]
m3<-m2[order(m2$subtype),]
sam<-data.frame(row.names = m3$sample,subtype=m3$subtype)
ha <- HeatmapAnnotation(df=sam,show_annotation_name = TRUE)
rownames(m3)<-m3$sample
mplot<-as.data.frame(t(m3))
mplot<-mplot[-c(1,11),]
mplot[mplot==0]<-NA
# ==================== Fig. 3b ====================
# Nine-gene oncoplot with subtype annotation; mutations are shown in black.
col_order<-colnames(mplot)
row_order<-c('STK11','KEAP1','TPTE','KCNU1','SPEG','SLC39A12','GRIK2','MMP16','AFF2')
oncoPrint(mplot,alter_fun = alter_fun,bottom_annotation = ha,column_order=col_order,row_order = row_order)


# ==============================================================================
# Part 5. Compare response-associated marker expression (Fig. 3c)
# ==============================================================================
# Average log2(TPM + 1) across markers using a fixed divisor of 17.
gep<-read.table('immunotherapy_biomarkers.txt')
load('exp_tpm.RData')
rownames(exp_tpm)<-exp_tpm[,2]
exp_tpm<-exp_tpm[,-c(1,2)]
exp_tpm<-log2(exp_tpm+1)
kras_sample<-read.csv('kras_sample.csv')
tpmkras<-subset(exp_tpm,select = kras_sample$sample)
save(tpmkras,file='tpmkras.RData')

gep_tpm<-tpmkras[which(rownames(tpmkras) %in% gep$V1),]
save(gep_tpm,file='gep_tpm.RData')

cold_id<-kras_sample[which(kras_sample$subtype=='cold'),]
gep_cold<-subset(gep_tpm,select=cold_id$sample)
cold<-colSums(gep_cold)/17

hot_id<-kras_sample[which(kras_sample$subtype=='hot'),]
gep_hot<-subset(gep_tpm,select=hot_id$sample)
hot<-colSums(gep_hot)/17

data<-data.frame(subtype=c(rep('hot',72),rep('cold',66)),gep=c(hot,cold))
data$subtype<-factor(data$subtype,levels = c('hot','cold'))
colnames(data)[2]<-'aver_exp'
immunotherapydata<-data
save(immunotherapydata,file='immunotherapydata.RData')
load(file='immunotherapydata.RData')
write.csv(immunotherapydata,file='aver_exp_of_immunotherapy_biomarker.csv')
# ==================== Fig. 3c ====================
# Density of mean log2(TPM + 1); displayed P value is fixed.
library(ggplot2)
library(ggpubr)
p<-ggplot(immunotherapydata,aes(x=aver_exp))
p+geom_density(aes(fill=subtype),alpha=0.4)+labs(x='The average expression of markers indicating immunotherapy response',y='density')+theme_classic()+
  annotate("text",x=7,y=0.5,label='p<2.2e-16')
# Test one marker's normality and compare subtype mean scores.
shapiro.test(as.matrix(gep_tpm[4,]))
t.test(aver_exp~subtype,data)

# ==============================================================================
# Part 6. GO biological-process enrichment: cold versus hot
# ==============================================================================
# Convert symbols to Entrez IDs and enrich GO biological processes.
library(clusterProfiler)
library(org.Hs.eg.db)
gobp<-function(deg){
  deg_entr<-bitr(geneID = deg,
                 fromType = 'SYMBOL',
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
  gombp<- enrichGO(gene = deg_entr$ENTREZID,
                   keyType = "ENTREZID",
                   OrgDb=org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500,
                   qvalueCutoff = 0.2,
                   readable = TRUE)
  gores<-gombp@result
  return(gores)
}

ch<-read.csv(file= "all_deseq2_results_cold_vs_hot.csv")
ch_up<-ch[which(ch$log2FoldChange>log2(1.5) & ch$padj<0.05),]
ch_down<-ch[which(ch$padj<0.05 & ch$log2FoldChange< - log2(1.5)),]
# Enrich genes increased and decreased in cold tumors.
gobp_chup<-gobp(ch_up$gene_name)
write.csv(gobp_chup,file='GObp_of_ch_up_degs.csv')

gobp_chdown<-gobp(ch_down$gene_name)
write.csv(gobp_chdown,file='GObp_of_ch_down_degs.csv')

# Combine the first 10 terms from each direction.
gobp_chupres1<-gobp_chup[1:10,]
gobp_chdownres1<-gobp_chdown[c(1:10),]
gobp_chres<-rbind(gobp_chdownres1,gobp_chupres1)

# Despite its name, logqvalue uses raw P values.
gobp_chres$logqvalue<-c(-log10(gobp_chdownres1$pvalue),-log10(gobp_chupres1$pvalue))
gobp_chres$group<-factor(c(rep('down-regulation',10),rep('up-regulation',10)),levels = c('up-regulation','down-regulation'))
colnames(gobp_chres)[11]<-'cold vs hot'
colnames(gobp_chres)[2]<-'Enrichment of Gene Ontology Biological Process'
write.csv(gobp_chres,file='GObp_of_ch_degs.csv')

library(ggpubr)
library(ggplot2)
save(gobp_chres,file='gobp_chres.RData')
load(file='gobp_chres.RData')
# ==================== Fig. 3e ====================
ggbarplot(gobp_chres,x='Enrichment of Gene Ontology Biological Process',y='logqvalue',fill='cold vs hot',color='white',rotate = T,sort.val = "asc",ggtheme = theme_minimal())

# ==============================================================================
# Part 7. TIMER2 immune-infiltration comparison (Supplementary Fig. 1a)
# ==============================================================================
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

timer2_kras<-read.csv(paste0(WORK_DIR,'/Results/TIMER2_KRAS_mutation_table.csv'))
timer2_kras_luad<-timer2_kras %>%
  filter(cancer=='LUAD (n=515)')
colnames(timer2_kras_luad)<-c("cancer","infiltrates","p","log2FC","x")

timer2_kras_luad <- timer2_kras_luad %>%
  mutate(
    cell_type_main = str_remove(infiltrates, "_.*$"),
    method = str_remove(infiltrates, "^[^_]+_"),
    log2FC = ifelse(is.infinite(log2FC), NA, log2FC)
  ) %>%
  group_by(cell_type_main) %>%
  mutate(
    n_method_p005 = n_distinct(method[p < 0.05]),
    n_method_p010 = n_distinct(method[p < 0.10]),
    n_tool_valid = n_distinct(method[!is.na(log2FC)])
  ) %>%
  ungroup()


data<-timer2_kras_luad %>%
  filter(n_tool_valid>1)
data<-timer2_kras_luad 
# NOTE: This assignment overwrites the filter above.
# Order cell types by mean log2FC.
cell_order <- data %>%
  group_by(cell_type_main) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_log2FC) %>%
  pull(cell_type_main)

plot_df <- data %>%
  select(cell_type_main, method, log2FC, p) %>%
  distinct() %>%
  complete(cell_type_main, method) %>%
  mutate(
    cell_type_main = factor(cell_type_main, levels = cell_order),
    sig_label = ifelse(!is.na(p) & p < 0.05, "*", "")
  )

# ==================== Supplementary Fig. 1a ====================
# TIMER2 heatmap: fill = log2FC; * = P < 0.05; gray = missing value.
p <- ggplot(plot_df, aes(x = cell_type_main, y = method, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 4) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey80",
    name = "log2FC"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p

