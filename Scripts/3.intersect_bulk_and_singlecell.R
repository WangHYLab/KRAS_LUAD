# Figure overview:
#   Fig. 4b: Bulk RNA-seq expression heatmap
#   Fig. 4c: Single-cell expression heatmap
#   Fig. 4d: KEGG enrichment analysis
#   Fig. 4e: Immune-suppressive CIIR genes after 24-h KRAS inhibition or 72-h KRAS inhibition
#   Fig. 4f: Immune-promoting CIIR genes after 24-h KRAS inhibition or 72-h KRAS inhibition
#   Fig. 4g: 
setwd('path/to/data')
# ==============================================================================
# Part 1. Identify shared immune-associated genes
# ==============================================================================

# Bulk RNA-seq: cold versus hot tumors.
chdeseq2<-read.csv('all_deseq2_results_cold_vs_hot.csv')# Cold versus hot.
chdeg_up<-chdeseq2[which(chdeseq2$log2FoldChange>log2(1.5) & chdeseq2$padj<0.05),]
chdeg_down<-chdeseq2[which(chdeseq2$log2FoldChange< -log2(1.5) & chdeseq2$padj<0.05),]

# scRNA-seq: genes correlated with immune-cell abundance.
sc<-read.csv('sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv')
sc_inhibit<-sc[which(sc$r<0 & sc$p<0.05),]
sc_activate<-sc[which(sc$r>0 & sc$p<0.05),]

# Intersect bulk DEGs with scRNA-seq correlation results.
A<-intersect(chdeg_down$gene_name,sc_inhibit$X)
immune_inhibit<-intersect(chdeg_up$gene_name,sc_inhibit$X)
immune_activate<-intersect(chdeg_down$gene_name,sc_activate$X)
ch_data<-chdeseq2[which(chdeseq2$gene_name %in% c(immune_inhibit,immune_activate)),]
colnames(ch_data)<-c('gene_name','baseMean_ch','log2FoldChange_ch','lfcSE_ch','stat_ch','pvalue_ch','padj_ch')
sc_data<-sc[which(sc$X %in% c(immune_inhibit,immune_activate)),-1]
colnames(sc_data)<-c('r_sc','p_sc','immune_percentage_sc','malignant_cell_exp_mean_sc','gene_name','p.adj_sc')
immune_genes_data<-merge(ch_data,sc_data,by='gene_name')
write.csv(immune_genes_data,file = 'immune_genes_data.csv',row.names = F)
immune_genes_data<-read.csv(file = 'immune_genes_data.csv')
# ==============================================================================
# Part 2. Annotate shared genes
# ==============================================================================

# Retrieve gene descriptions from Ensembl BioMart.
library(biomaRt)
my_mart_h <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info_h <- getBM(attributes=c("external_gene_name","description"),filters = "external_gene_name",
                     values=immune_genes_data$gene_name,
                     mart = my_mart_h)
gene_info_h1<-gene_info_h[order(gene_info_h$external_gene_name),]
colnames(gene_info_h1)<-c('gene_name','gene_description')
o1<-merge(gene_info_h1,immune_genes_data,by='gene_name',all=T)
o1$immune<-rep('negative',nrow(o1))
o1[which(o1$r_sc>0),'immune']<-'positive'
write.csv(o1,file='immune_genes_basical_info.csv',row.names=F)
o1<-read.csv(file='immune_genes_basical_info.csv')
o1$nvalid_sample <- sapply(o1$malignant_cell_exp_mean_sc, function(x) {
  sum(!is.na(as.numeric(strsplit(x, ",")[[1]])))
})
write.csv(o1,file='overlap_65_genes_of_bulk_and_sc.csv',row.names=F)



# ==============================================================================
# Part 3. Visualize expression in bulk and scRNA-seq data
# ==============================================================================

# Bulk RNA-seq heatmap of shared genes.
immune_genes<-read.csv('immune_genes_basical_info.csv')
library(DESeq2)
library(pheatmap)
load(file='xx/Rdata/dds2_ch.RData')
ch_dds<-counts(dds2_ch,normalized=T)
ch_exp<-ch_dds[which(rownames(ch_dds) %in% immune_genes$gene_name),]
ch_exp<-log2(ch_exp+1)
save(ch_exp,file='ch_exp.RData')
kras_sample<-read.csv('kras_sample.csv')
rownames(kras_sample)<-kras_sample$sample
annotation_col_ch<-kras_sample[,-1,drop=F]

gene<-unique(immune_genes[,c('gene_name','immune')])
#gene1<-gene[which(gene$gene_name %in% rownames(ch_exp)),]
rownames(gene)<-gene$gene_name
annotation_row_ch<-gene[,-1,drop=F]
annotation_row_ch<-annotation_row_ch[order(annotation_row_ch$immune,decreasing = T),,drop=F]

count1_ch<-ch_exp[rownames(annotation_row_ch),]
z_score <- function(x) {
  return((x - mean(x)) / sd(x))
}
#count_ch <- as.data.frame(t(scale(t(count1_ch)))) # Equivalent to the custom z_score function.

count_ch <- as.data.frame(t(apply(count1_ch, 1, z_score)))
#count_ch1<-count_ch[which(rownames(count_ch) %in% data1$gene_name),]
colors <- colorRampPalette(c("blue", "white", "red"))(100)
max_val <- max(abs(count_ch))  # Maximum absolute matrix value.
breaks <- seq(-max_val, max_val, length.out = 101) 
count_ch_s<-count_ch[rownames(data5)[-1],]
###########Fig4b##########
p_ch<-pheatmap(count_ch_s,color = colors, breaks = breaks,show_colnames = F,show_rownames = T,cluster_rows = F,cluster_cols = F,fontsize_col=5,fontsize_row=9,annotation_col = annotation_col_ch,annotation_row = annotation_row_ch,gaps_row = c(12),gaps_col = c(66))
# scRNA-seq heatmap of shared genes.
library(pheatmap)
sc<-read.csv('sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csv')
immune_genes<-read.csv('immune_genes_basical_info.csv')
sc_immune_genes<-sc[which(sc$X %in% immune_genes$gene_name),]
df<-sc_immune_genes
library(tidyverse)
n_samples <- max(str_count(df$immune_percentage, ",") + 1)  # Number of samples.
sample_names <- paste0("Sample_", 1:n_samples)

df_split <- df %>%
  mutate(across(malignant_cell_exp_mean, ~ str_split(.x, ","))) %>%
  unnest_wider(malignant_cell_exp_mean, names_sep = "_") 

immune_values <- str_split(df[1, 4], ",")[[1]]  # Split immune-cell percentages into a vector.
immune_percentage <- c("Immune_Percentage", immune_values)
data<-rbind(df_split[,c(27,5:26)],immune_percentage)
colnames(data)<-c("gene_name",sample_names)
#data[data == "NA"] <- NA  # Convert the string "NA" to missing values.
#data_na<- data[rowSums(is.na(data)) <= 11, ]
data1<-data
#data1<-data_na
# Create row annotations for gene categories.
annotation_row <- annotation_row_ch[which(rownames(annotation_row_ch) %in% data1$gene_name),,drop=F] # Match annotation and matrix row names.
annotation_row['Immune_Percentage',]<-'Immune_Percentage'
annotation_row$immune <- factor(annotation_row$immune, levels = c("Immune_Percentage", "positive", "negative"))
annotation_row<-annotation_row[order(annotation_row$immune),,drop=F]

ann_colors <- list(
  Category = c(Immune_Percentage = "purple", positive = "red", negative = "blue")
)
data2<-as.data.frame(data1[,-1])
data3 <- as.data.frame(lapply(data2, function(x) as.numeric(as.character(x))))
rownames(data3)<-data1$gene_name
order_imm<-as.data.frame(t(data3['Immune_Percentage',]))
order_imm<-order_imm[order(order_imm$Immune_Percentage),,drop=F]

data3_sorted <- data3[,rownames(order_imm)]  # Reorder samples by immune-cell percentage.
data3_sorted1<-data3_sorted[rownames(annotation_row),]
data3_sorted1$immune<-annotation_row$immune
na_count <- apply(data3_sorted1, 1, function(x) sum(is.na(x)))
data4 <- data3_sorted1[order(na_count), ]
data5<-data4[order(data4$immune),]

# ###########Fig4c##########
p_sc<-pheatmap(
  data5[,-23],
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_row = annotation_row,  # Separate immune percentage, positive genes, and negative genes.
  annotation_colors = ann_colors,
  na_col = "gray",
  gaps_row = c(1,13),
  fontsize_row = 9,
  show_colnames = FALSE
)
grid.arrange(p_ch$gtable, p_sc$gtable, ncol = 2, widths = c(4.5,5))  # Arrange the heatmaps horizontally.



# ==============================================================================
# Part 4. Perform KEGG enrichment analysis
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Run KEGG enrichment for a gene set.
kegg<-function(deg){
  deg_entr<-bitr(geneID = deg,
                 fromType = 'SYMBOL',
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
  kegge<- enrichKEGG(gene = deg_entr$ENTREZID,  # Converted Entrez IDs.
                     keyType = "kegg",
                     organism = "hsa",
                     pvalueCutoff = 0.05,  # P-value cutoff.
                     pAdjustMethod = "BH",  # Multiple-testing correction.
                     minGSSize = 10,   # Minimum gene-set size.
                     maxGSSize = 500,  # Maximum gene-set size.
                     qvalueCutoff = 0.2) 
  keggres<-kegge@result
  keggres$geneID <- sapply(keggres$geneID, function(ids) {
    genes <- unlist(strsplit(ids, "/"))
    name_df <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    symbol_vec <- name_df$SYMBOL[match(genes, name_df$ENTREZID)]
    paste(symbol_vec, collapse = "/")
  })
  return(keggres)
}
# Load annotated immune-associated genes.
immune_genes<-read.csv('immune_genes_basical_info.csv')

# Enrich positive and negative gene sets separately.
kegg_pos<-kegg(deg = immune_genes[which(immune_genes$immune=='positive'),]$gene_name)
kegg_neg<-kegg(deg =immune_genes[which(immune_genes$immune=='negative'),]$gene_name)
write.csv(kegg_pos,file='immune_genes_kegg_pos.csv',row.names = F)
write.csv(kegg_neg,file='immune_genes_kegg_neg.csv',row.names = F)
# Select pathways for visualization.
kegg_pos$Group<-rep('positive',nrow(kegg_pos))
kegg_neg$Group<-rep('negative',nrow(kegg_neg))
kegg<-rbind(kegg_neg[1:5,],kegg_pos[1:3,])
###########Fig4d##########
ggplot(kegg, aes(x = GeneRatio, 
                    y = reorder(Description, GeneRatio), 
                    size = Count, 
                    color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(Group ~ ., scales = "free", space = "free") +  # Stack groups vertically.
  labs(title = "KEGG Enrichment of immune-related genes",
       x = "Gene Ratio", y = NULL, color = "Adjust P-value") +
  theme_bw() +
  theme(strip.text.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust=0.5, size=12),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7))




# ==============================================================================
# Part 5. examine CIIR genes in KRAS-inhibitor perturbation datasets
# ==============================================================================

# H358 cells treated with ARS-1620 for 24 h or adagrasib (MRTX849) for 72 h.
WORK_DIR<-"path/to/KRAS_LUAD/"
library(tidyverse)
library(readr)
library(purrr)

library(patchwork)
# Extract upregulated and downregulated genes.
get_deg <- function(df,
                    gene_col,
                    lfc_col = "log2FoldChange",
                    padj_col = "padj",
                    lfc_cut = log2(1.5),
                    padj_cut = 0.05) {
  
  df <- df[!is.na(df[[padj_col]]), ]
  
  list(
    up   = unique(df[df[[lfc_col]] >  lfc_cut & df[[padj_col]] < padj_cut, gene_col]),
    down = unique(df[df[[lfc_col]] < -lfc_cut & df[[padj_col]] < padj_cut, gene_col])
  )
}



# Build the GSE275534 feature-count matrix.
data_dir <- "path/to/GSE275534_RAW/"

files <- list.files(
  data_dir,
  pattern = "featurecounts_genes.tsv.gz$",
  full.names = TRUE
)

length(files)
sample_names <- basename(files) %>%
  sub("\\.featurecounts_genes.tsv.gz$", "", .)

sample_names
read_fc <- function(file) {
  df <- read_tsv(
    file,
    comment = "#",
    col_types = cols()
  )
  
  df %>%
    dplyr::select(
      gene_id = Geneid,
      count = ncol(df)
    )
}
count_list <- lapply(files, read_fc)

count_mat <- reduce(
  count_list,
  full_join,
  by = "gene_id"
)

colnames(count_mat)[-1] <- sample_names
count_mat<-as.data.frame(count_mat)
save(count_mat,file='./Rdata/count_mat.RData')
load(file='./Rdata/count_mat.RData')
rownames(count_mat)<-count_mat$gene_id
count_mat<-count_mat[,-1]

# Compare 72-h adagrasib treatment with control.
data_72h<-count_mat[,c(1:4,9:12)]

sample_info_72h <- data.frame(
  row.names = colnames(data_72h),
  group = factor(c("control","control","control","control","72h","72h","72h","72h"), levels = c("control", "72h"))
)
library(DESeq2)

dds_72h <- DESeqDataSetFromMatrix(countData =data_72h,
                                  colData = sample_info_72h,
                                  design = ~ group)
dds_72h <- dds_72h[rowSums(counts(dds_72h)) > 10, ]
dds_72h <- DESeq(dds_72h)
res_72h <- results(dds_72h, contrast = c("group","72h","control"))
res_72h <- res_72h[order(res_72h$pvalue), ]
res_df_72h <- as.data.frame(res_72h)
res_df_72h$gene<-rownames(res_df_72h)

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- res_df_72h$gene

anno <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)
res_df_72h1<-merge(res_df_72h,anno,by.x='gene',by.y='ensembl_gene_id')
write.csv(res_df_72h1,file='./Results/GSE275534_h358_adag_72h_vs_control_deseq2.csv',row.names = F)


# Prepare GSE275534 differential-expression gene sets.
GSE275534_h358_72h<-read.csv(paste0(WORK_DIR,"Results/GSE275534_h358_adag_72h_vs_control_deseq2.csv"))
GSE275534_h358_72h<-GSE275534_h358_72h[which(GSE275534_h358_72h$hgnc_symbol!=''),]
GSE275534_h358_72h_deg<-get_deg(GSE275534_h358_72h,gene_col = 'hgnc_symbol',lfc_cut=0)
# Load immune-promoting and immune-suppressive CIIR genes.
ciir_genes<-read.csv(paste0(WORK_DIR,"Results/overlap_65_genes_of_bulk_and_sc.csv"))
immune_negative<-ciir_genes[which(ciir_genes$immune=='negative'),]$gene_name
immune_positive<-ciir_genes[which(ciir_genes$immune=='positive'),]$gene_name


GSE275534_h358_72h_gene_sets <- list(
  KRAS_inhibitor_up = GSE275534_h358_72h_deg$up,
  KRAS_inhibitor_down = GSE275534_h358_72h_deg$down,
  immune_positive   = immune_positive,
  immune_negative   = immune_negative
)

library(readr)
# Compare 24-h ARS-1620 treatment with control in GSE183516.
GSE183516<-read_tsv(file=paste0("xx/Data/GSE183516_raw_counts_GRCh38/GSE183516_raw_counts_GRCh38.p13_NCBI.tsv"))
GSE183516_h358<-GSE183516[,c('GeneID','GSM5558640','GSM5558641','GSM5558642','GSM5558649','GSM5558650','GSM5558651')]

gene_info<-read_tsv(file=paste0("x/Data/Human.GRCh38.p13.annot.tsv"))
data_24h <- merge(GSE183516_h358, gene_info[,c(1,2)],
             by = "GeneID")
data_24h<-data_24h[!duplicated(data_24h$Symbol),]
rownames(data_24h)<-data_24h$Symbol
data_24h<-data_24h[,-c(1,8)]
sample_info_24h <- data.frame(
  row.names = colnames(data_24h),
  group = factor(c("control","control","control","24h","24h","24h"), levels = c("control", "24h"))
)
library(DESeq2)

dds_24h <- DESeqDataSetFromMatrix(countData =data_24h,
                                  colData = sample_info_24h,
                                  design = ~ group)
dds_24h <- dds_24h[rowSums(counts(dds_24h)) > 10, ]
dds_24h <- DESeq(dds_24h)
res_24h <- results(dds_24h, contrast = c("group","24h","control"))
res_24h <- res_24h[order(res_24h$pvalue), ]
res_df_24h <- as.data.frame(res_24h)
res_df_24h$gene<-rownames(res_df_24h)
write.csv(res_df_24h,paste0(WORK_DIR,"Results/GSE183516_h358_ars1620_24h_vs_control_deseq2.csv"))

GSE183516_h358_24h_deg<-get_deg(res_df_24h,gene_col = 'gene',lfc_cut=0)
GSE183516_h358_24h_gene_sets <- list(
  KRAS_inhibitor_up = GSE183516_h358_24h_deg$up,
  KRAS_inhibitor_down = GSE183516_h358_24h_deg$down,
  immune_positive   = immune_positive,
  immune_negative   = immune_negative
)



# Select significant CIIR genes and compare their log2 fold changes.
GSE275534_h358_72h<-read.csv(paste0(WORK_DIR,"Results/GSE275534_h358_adag_72h_vs_control_deseq2.csv"))
GSE275534_h358_72h<-GSE275534_h358_72h[!duplicated(GSE275534_h358_72h$hgnc_symbol),]
rownames(GSE275534_h358_72h)<-GSE275534_h358_72h$hgnc_symbol
GSE275534_h358_72h<-GSE275534_h358_72h[which(GSE275534_h358_72h$hgnc_symbol!=''),]

GSE275534_h358_72h_deg<-GSE275534_h358_72h[which(GSE275534_h358_72h$padj<0.05),]



GSE183516_h358_24h<-read.csv(paste0(WORK_DIR,"Results/GSE183516_h358_ars1620_24h_vs_control_deseq2.csv"),row.names = 1)
GSE183516_h358_24h<-GSE183516_h358_24h[which(GSE183516_h358_24h$gene!=''),]
GSE183516_h358_24h_deg<-GSE183516_h358_24h[which(GSE183516_h358_24h$padj<0.05),]


# Reload CIIR gene groups.
ciir_genes<-read.csv(paste0(WORK_DIR,"Results/fig4.overlap_65_genes_of_bulk_and_sc.csv"))
immune_negative<-ciir_genes[which(ciir_genes$immune=='negative'),]$gene_name
immune_positive<-ciir_genes[which(ciir_genes$immune=='positive'),]$gene_name

# Identify overlapping genes for each treatment and immune group.
h358_72h_negative<-GSE275534_h358_72h_deg[which(GSE275534_h358_72h_deg$hgnc_symbol %in% immune_negative),]
h358_72h_negative$immune_group<-rep('suppressive',nrow(h358_72h_negative))
h358_72h_positive<-GSE275534_h358_72h_deg[which(GSE275534_h358_72h_deg$hgnc_symbol %in% immune_positive),]
h358_72h_positive$immune_group<-rep('promoting',nrow(h358_72h_positive))

h358_72h_genes<-rbind(h358_72h_negative,h358_72h_positive)
colnames(h358_72h_genes)<-c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","gene","immune_group")

h358_24h_negative<-GSE183516_h358_24h_deg[which(GSE183516_h358_24h_deg$gene %in% immune_negative),]
h358_24h_negative$immune_group<-rep('suppressive',nrow(h358_24h_negative))
h358_24h_positive<-GSE183516_h358_24h_deg[which(GSE183516_h358_24h_deg$gene %in% immune_positive),]
h358_24h_positive$immune_group<-rep('promoting',nrow(h358_24h_positive))

h358_24h_genes<-rbind(h358_24h_negative,h358_24h_positive)


# Plot treatment-induced log2 fold changes for one immune group.
library(dplyr)
library(ggplot2)

plot_log2FC_single <- function(data,
                               group_to_plot = c("promoting", "suppressive"),
                               group_label = NULL,
                               fill_color = "grey35",
                               border_color = NA) {
  
  group_to_plot <- match.arg(group_to_plot)
  
  data_plot <- data %>%
    filter(immune_group == group_to_plot) %>%
    arrange(desc(log2FoldChange)) %>%
    mutate(gene = factor(gene, levels = gene))
  
  if (is.null(group_label)) {
    group_label <- c(
      "promoting" = "immune-promoting",
      "suppressive" = "immune-suppressive"
    )[group_to_plot]
  }
  
  p <- ggplot(data_plot, aes(x = gene, y = log2FoldChange)) +
    geom_col(width = 0.7, fill = fill_color, color = border_color) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(
      title = "",
      x = group_label,
      y = "log2FoldChange (KRAS inhibitor vs control)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 10),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6)
    )
  
  return(p)
}
########Fig4e-f#####
p1 <- plot_log2FC_single(
  h358_24h_genes,
  group_to_plot = "promoting",
  group_label = "immune-promoting CIIR genes",
  fill_color = "grey80",
  border_color = "black"
)
p2 <- plot_log2FC_single(
  h358_24h_genes,
  group_to_plot = "suppressive",
  group_label = "immune-suppressive CIIR genes",
  #fill_color = "grey35"
  fill_color = "grey80",
  border_color = "black"
)
p3 <- plot_log2FC_single(
  h358_72h_genes,
  group_to_plot = "promoting",
  group_label = "immune-promoting CIIR genes",
  fill_color = "grey80",
  border_color = "black"
)
p4 <- plot_log2FC_single(
  h358_72h_genes,
  group_to_plot = "suppressive",
  group_label = "immune-suppressive CIIR genes",
  #fill_color = "grey35"
  fill_color = "grey80",
  border_color = "black"
)
grid.arrange(p2,p4,p1,p3,ncol=4,widths=c(4,2,2,1.5))

####################################################################
###################Fig4g and supplementary figure 3-4#########
# 1. gene groups
# Genes shared by the two studies
genes_suppressive_both <- c("PDE4D", "KCNF1", "DEPDC1")
genes_promoting_both <- c("CYP1B1", "MMP7", "CX3CL1", "PLAT")

# 2. Analysis and plotting helpers ----
# Expression: genes x samples. TIMER2: cell_type column followed by sample columns.
get_gene_vs_all_celltypes <- function(exp_mat, cell_mat, gene_name) {
  if (!(gene_name %in% rownames(exp_mat))) {
    stop("gene not found")
  }
  
  gene_vec <- as.numeric(exp_mat[gene_name, ])
  names(gene_vec) <- colnames(exp_mat)
  
  res_list <- lapply(seq_len(nrow(cell_mat)), function(i) {
    cell_name <- cell_mat$cell_type[i]
    cell_vec <- as.numeric(cell_mat[i, -1])
    names(cell_vec) <- colnames(cell_mat)[-1]
    
    common_samples <- intersect(names(gene_vec), names(cell_vec))
    gene_use <- gene_vec[common_samples]
    cell_use <- cell_vec[common_samples]
    
    keep <- complete.cases(gene_use, cell_use)
    gene_use <- gene_use[keep]
    cell_use <- cell_use[keep]
    
    # Skip tests with fewer than three complete pairs or no variation.
    if (length(gene_use) < 3 || length(unique(gene_use)) < 2 || length(unique(cell_use)) < 2) {
      return(data.frame(
        gene = gene_name,
        cell_type = cell_name,
        cor = NA,
        pvalue = NA,
        direction = NA,
        n = length(gene_use)
      ))
    }
    
    cor_res <- cor.test(gene_use, cell_use, method = "spearman")
    cor_val <- unname(cor_res$estimate)
    
    data.frame(
      gene = gene_name,
      cell_type = cell_name,
      cor = cor_val,
      pvalue = cor_res$p.value,
      direction = ifelse(cor_val > 0, "positive", "negative"),
      n = length(gene_use)
    )
  })
  
  # Apply BH correction across all TIMER2 indicators within each gene.
  res_df <- bind_rows(res_list) %>%
    mutate(
      padj = p.adjust(pvalue, method = "BH")
    ) %>%
    arrange(pvalue)
  
  return(res_df)
}

# Heatmaps retain cell types with valid correlations from at least two methods.
get_timer2_gene_analysis <- function(gene_name, exp_mat, cell_mat,
                                    direction_name) {
  if (!(gene_name %in% rownames(exp_mat))) {
    warning(sprintf("Skipping %s: gene not found in kras_exp.", gene_name))
    return(NULL)
  }

  res <- get_gene_vs_all_celltypes(
    exp_mat = exp_mat,
    cell_mat = cell_mat,
    gene_name = gene_name
  )

  res_plot <- res %>%
    mutate(
      cell_type_main = str_remove(cell_type, "_[^_]+$"),
      method = str_extract(cell_type, "[^_]+$"),
      cor = ifelse(is.infinite(cor), NA, cor)
    ) %>%
    group_by(cell_type_main) %>%
    mutate(
      n_tool_valid = n_distinct(method[!is.na(cor)])
    ) %>%
    ungroup()

  res_plot_filter <- res_plot %>%
    filter(n_tool_valid > 1)

  if (nrow(res_plot_filter) == 0) {
    warning(sprintf("No cell types with at least two valid methods for %s.", gene_name))
    return(list(res = res, plot_df = NULL, p = NULL))
  }

  # Order by decreasing number of negative methods, then increasing mean correlation.
  cell_order <- res_plot_filter %>%
    group_by(cell_type_main) %>%
    summarise(
      n_neg_method = n_distinct(method[!is.na(cor) & cor < 0]),
      mean_cor = mean(cor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_neg_method), mean_cor) %>%
    pull(cell_type_main)

  # Heatmap stars use raw P < 0.05; bubble sizes use BH-adjusted P.
  plot_df <- res_plot_filter %>%
    dplyr::select(cell_type_main, method, cor, pvalue, padj) %>%
    distinct() %>%
    complete(cell_type_main, method) %>%
    mutate(
      sig_label = case_when(
        !is.na(pvalue) & pvalue < 0.05 ~ "*",
        TRUE ~ ""
      ),
      cell_type_main = factor(cell_type_main, levels = cell_order)
    )

  p <- ggplot(plot_df, aes(x = cell_type_main, y = method, fill = cor)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 4) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, na.value = "grey80", name = "Spearman\ncorrelation"
    ) +
    labs(
      title = gene_name,
      subtitle = paste(direction_name, "Both studies", sep = " | "),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  list(res = res, plot_df = plot_df, p = p)
}

# Arrange the gene heatmaps within each group.
arrange_timer2_group <- function(plots, ncol = 2, title = NULL) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) {
    message("No plots to arrange: ", title)
    return(NULL)
  }
  gridExtra::grid.arrange(
    grobs = plots,
    ncol = min(ncol, length(plots)),
    top = title
  )
}

# Match linked IDs as sets (A;B = B;A), keeping each RNA column as one sample.
normalize_cptac_sample_ids <- function(ids) {
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop("Missing CPTAC sample identifiers.")
  }
  vapply(strsplit(ids, ";", fixed = TRUE), function(x) {
    paste(sort(unique(trimws(x))), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}
timer2_gene_groups <- list(
  suppressive = list(both = genes_suppressive_both),
  promoting = list(both = genes_promoting_both)
)
# Shared cohort workflow; results are indexed by $direction$both$gene.
analyse_timer2_both <- function(exp_mat, cell_mat, cohort = NULL) {
  results <- lapply(names(timer2_gene_groups), function(direction_name) {
    genes <- timer2_gene_groups[[direction_name]]$both
    direction_label <- if (is.null(cohort)) direction_name else
      paste(cohort, direction_name, sep = " | ")
    gene_results <- lapply(genes, function(gene_name) {
      get_timer2_gene_analysis(gene_name, exp_mat, cell_mat, direction_label)
    })
    list(both = setNames(gene_results, genes))
  })
  setNames(results, names(timer2_gene_groups))
}

cptac_analysis_genes <- unique(c(genes_suppressive_both, genes_promoting_both))

# 3. TCGA inputs and shared-gene heatmaps ----
load(paste0(WORK_DIR,'Rdata/kras_tpm138.RData'))
kras_tpm<-kras_tpm138[!duplicated(kras_tpm138$gene_name),]
kras_tpm_coding<-kras_tpm[which(kras_tpm$genetype=='protein_coding'),]
rownames(kras_tpm_coding)<-kras_tpm_coding$gene_name
kras_exp<-kras_tpm_coding[,-c(1,2)]

matrix<-read.csv(paste0(WORK_DIR,'Results/TIMER2_estimation_matrix_KRAStpm.csv'))
colnames(matrix) <- gsub("\\.", "-", colnames(matrix))

timer2_results <- analyse_timer2_both(kras_exp, matrix)

p_list <- lapply(timer2_results, function(group) {
  list(both = lapply(group$both, function(result) result$p))
})

p_group <- list(suppressive = list(), promoting = list())

p_group$promoting$both <- arrange_timer2_group(p_list$promoting$both, ncol = 1, title = "both promoting")
p_group$suppressive$both <- arrange_timer2_group(p_list$suppressive$both, ncol = 1, title = "both suppressive")

# 4. Optional CPTAC TPM preparation ----
# GDC: https://portal.gdc.cancer.gov/analysis_page?app=CohortBuilder&tab=genomic_filters
# Select CPTAC-3, bronchus and lung, and adenomas and adenocarcinomas; download masked somatic mutations.
# Original selection: 236 cases; 113 matched the manifest, including 68 cases with nonsynonymous KRAS mutations (69 MAF files).

# load(file=paste0(WORK_DIR,'Rdata/exp_tpm_cptac.RData'))
# kras_cptac<-read.csv(paste0(WORK_DIR,'Tables/KRAS_mutant_RNA_samples_cptac.csv'),
#                      check.names = FALSE, stringsAsFactors = FALSE)
# 
# cptac_tpm_columns <- setdiff(colnames(exp_tpm), c("genetype", "gene_name"))
# cptac_tpm_keys <- normalize_cptac_sample_ids(cptac_tpm_columns)
# kras_cptac_keys <- normalize_cptac_sample_ids(kras_cptac$sample_submitter_ids)
# if (anyDuplicated(cptac_tpm_keys) || anyDuplicated(kras_cptac_keys) ||
#     anyDuplicated(kras_cptac$case_submitter_id)) {
#   stop("Duplicate CPTAC sample sets or cases; review before extracting TPM.")
# }
# cptac_match <- match(kras_cptac_keys, cptac_tpm_keys)
# if (anyNA(cptac_match)) {
#   stop(paste("CPTAC KRAS samples missing from TPM:",
#              paste(kras_cptac$sample_submitter_ids[is.na(cptac_match)], collapse = ", ")))
# }
# 
# # Preserve mutation-table order; each case has one RNA file and is labeled by case ID.
# kras_cptac_sample_map <- kras_cptac
# kras_cptac_sample_map$tpm_source_column <- cptac_tpm_columns[cptac_match]
# kras_cptac_sample_map$tpm_output_column <- kras_cptac$case_submitter_id
# stopifnot(
#   all(kras_cptac_sample_map$sample_types == "Primary Tumor"),
#   all(kras_cptac_sample_map$KRAS_status_linked_sample ==
#         "KRAS_nonsynonymous_in_linked_MAF"),
#   all(vapply(exp_tpm[kras_cptac_sample_map$tpm_source_column], is.numeric, logical(1)))
# )
# 
# cptac_gene_names <- trimws(as.character(exp_tpm$gene_name))
# cptac_named_gene <- !is.na(cptac_gene_names) & nzchar(cptac_gene_names)
# cptac_coding_gene <- cptac_named_gene & !is.na(exp_tpm$genetype) &
#   exp_tpm$genetype == "protein_coding"
# 
# # Keep the first protein-coding row per symbol, retaining untransformed TPM.
# cptac_keep_gene <- cptac_coding_gene
# cptac_keep_gene[cptac_coding_gene] <-
#   !duplicated(cptac_gene_names[cptac_coding_gene])
# cptac_selected_tpm <- as.matrix(
#   exp_tpm[cptac_keep_gene, kras_cptac_sample_map$tpm_source_column, drop = FALSE]
# )
# stopifnot(all(is.finite(cptac_selected_tpm)), all(cptac_selected_tpm >= 0))
# 
# kras_tpm_cptac <- cptac_selected_tpm
# rownames(kras_tpm_cptac) <- cptac_gene_names[cptac_keep_gene]
# colnames(kras_tpm_cptac) <- kras_cptac_sample_map$tpm_output_column
# stopifnot(
#   !anyDuplicated(rownames(kras_tpm_cptac)),
#   identical(colnames(kras_tpm_cptac), kras_cptac$case_submitter_id),
#   isTRUE(all.equal(unname(colSums(kras_tpm_cptac)),
#                    unname(colSums(cptac_selected_tpm))))
# )
# 
# save(kras_tpm_cptac, kras_cptac_sample_map,
#      file = file.path(WORK_DIR, "Rdata/kras_tpm_cptac_TIMER2.RData"))
# 
# # Export TIMER2 input: gene symbols followed by sample TPM columns.
# write.table(data.frame(gene = rownames(kras_tpm_cptac), kras_tpm_cptac,
#                        check.names = FALSE),
#             file.path(WORK_DIR, "Results/TIMER2_CPTAC_KRAS_mutant_TPM.tsv"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# 5. CPTAC inputs and shared-gene heatmaps ----
load(file = file.path(WORK_DIR, "Rdata/kras_tpm_cptac_TIMER2.RData"))
cptac_matrix<-read.csv(paste0(WORK_DIR,'Results/TIMER2_CPTAC_KRAS_LUAD_estimation_matrix.csv'),
                       check.names = FALSE, stringsAsFactors = FALSE)

# Preserve sample IDs and align TIMER2 columns to the expression matrix.
if (!identical(names(cptac_matrix)[1], "cell_type") ||
    anyDuplicated(names(cptac_matrix)[-1]) || anyDuplicated(cptac_matrix$cell_type)) {
  stop("Invalid or duplicated CPTAC TIMER2 sample/cell-type identifiers.")
}
cptac_missing_samples <- setdiff(colnames(kras_tpm_cptac), names(cptac_matrix)[-1])
cptac_extra_samples <- setdiff(names(cptac_matrix)[-1], colnames(kras_tpm_cptac))
if (length(cptac_missing_samples) || length(cptac_extra_samples)) {
  stop(paste("CPTAC TPM/TIMER2 sample mismatch. Missing:",
             paste(cptac_missing_samples, collapse = ", "), "Extra:",
             paste(cptac_extra_samples, collapse = ", ")))
}
cptac_matrix <- cptac_matrix[, c("cell_type", colnames(kras_tpm_cptac)), drop = FALSE]
stopifnot(
  identical(names(cptac_matrix)[-1], colnames(kras_tpm_cptac)),
  all(vapply(cptac_matrix[-1], is.numeric, logical(1))),
  all(is.finite(as.matrix(cptac_matrix[-1])))
)
cptac_missing_genes <- setdiff(cptac_analysis_genes, rownames(kras_tpm_cptac))
if (length(cptac_missing_genes)) {
  stop(paste("Genes missing from CPTAC TPM:", paste(cptac_missing_genes, collapse = ", ")))
}

timer2_results_cptac <- analyse_timer2_both(kras_tpm_cptac, cptac_matrix, cohort = "CPTAC")

p_list_cptac <- lapply(timer2_results_cptac, function(group) {
  list(both = lapply(group$both, function(result) result$p))
})
p_group_cptac <- list(suppressive = list(), promoting = list())
p_group_cptac$promoting$both <- arrange_timer2_group(p_list_cptac$promoting$both, ncol = 1, title = "CPTAC | both promoting")
p_group_cptac$suppressive$both <- arrange_timer2_group(p_list_cptac$suppressive$both, ncol = 1, title = "CPTAC | both suppressive")

# 6. xCell immune score analysis and bubble plots for all 65 genes ----
# Use source-table groups: negative -> suppressive; positive -> promoting.
immune65_gene_table <- read.csv(
  file.path(WORK_DIR, "Results/fig4.overlap_65_genes_of_bulk_and_sc.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(!anyDuplicated(immune65_gene_table$gene_name),
          all(immune65_gene_table$immune %in% c("negative", "positive")))
immune65_gene_table <- immune65_gene_table %>%
  transmute(gene = gene_name,
            immune_group = factor(ifelse(immune == "negative", "suppressive", "promoting"),
                                  levels = c("promoting", "suppressive"))) %>%
  arrange(immune_group)
immune65_gene_order <- immune65_gene_table$gene
immune65_cohorts <- list(
  TCGA = list(exp = kras_exp, cells = matrix),
  CPTAC = list(exp = kras_tpm_cptac, cells = cptac_matrix)
)
immune65_results <- bind_rows(lapply(names(immune65_cohorts), function(cohort_name) {
  cohort_data <- immune65_cohorts[[cohort_name]]
  missing_genes <- setdiff(immune65_gene_order, rownames(cohort_data$exp))
  if (length(missing_genes)) {
    stop(paste(cohort_name, "missing genes:", paste(missing_genes, collapse = ", ")))
  }
  stopifnot(sum(cohort_data$cells$cell_type == "immune score_XCELL") == 1L,
            setequal(colnames(cohort_data$exp), names(cohort_data$cells)[-1]))
  # Test and adjust across all TIMER2 indicators before extracting xCell; ignore input-table padj.
  withCallingHandlers(
    bind_rows(lapply(immune65_gene_order, function(gene_name) {
      get_gene_vs_all_celltypes(cohort_data$exp, cohort_data$cells, gene_name) %>%
        filter(cell_type == "immune score_XCELL") %>%
        transmute(gene, cohort = cohort_name, cell_type, rho = cor, pvalue, padj, n)
    })),
    warning = function(w) {
      if (grepl("Cannot compute exact p-value with ties", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
})) %>%
  left_join(immune65_gene_table, by = "gene") %>%
  mutate(gene = factor(gene, levels = immune65_gene_order),
         cohort = factor(cohort, levels = c("TCGA", "CPTAC")),
         neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))
stopifnot(nrow(immune65_results) == 2L * nrow(immune65_gene_table),
          !anyDuplicated(immune65_results[c("gene", "cohort")]))

# Use common axis, color, and size scales across cohorts.
immune65_rho_limit <- max(0.1, ceiling(max(abs(immune65_results$rho), na.rm = TRUE) * 10) / 10)
immune65_size_limit <- max(immune65_results$neg_log10_padj, na.rm = TRUE)
p_immune65_base <- ggplot(immune65_results, aes(x = gene, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = rho, size = neg_log10_padj),
             shape = 21, color = "grey35", stroke = 0.35, na.rm = TRUE) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  scale_size_continuous(range = c(1.6, 7.5), limits = c(0, immune65_size_limit),
                        name = expression(-log[10](P[adj]))) +
  scale_y_continuous(limits = c(-immune65_rho_limit, immune65_rho_limit),
                     labels = scales::label_number(accuracy = 0.1),
                     expand = expansion(mult = 0.08)) +
  facet_grid(. ~ immune_group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Spearman rho") +
  guides(size = guide_legend(override.aes = list(fill = "white"))) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "right", plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8, color = "grey40"))

# Default order: increasing adjusted P within each cohort and group.
immune65_tcga_plot_data <- immune65_results %>%
  filter(cohort == "TCGA") %>%
  arrange(immune_group, padj, as.character(gene)) %>%
  mutate(gene = factor(as.character(gene), levels = as.character(gene)))
immune65_cptac_plot_data <- immune65_results %>%
  filter(cohort == "CPTAC") %>%
  arrange(immune_group, padj, as.character(gene)) %>%
  mutate(gene = factor(as.character(gene), levels = as.character(gene)))

p_tcga <- p_immune65_base %+% immune65_tcga_plot_data +
  labs(title = "Expression vs immune score_xCell TCGA (n = 138)")
p_cpatc <- p_immune65_base %+% immune65_cptac_plot_data +
  labs(title = "Expression vs immune score_xCell CPTAC (n = 68)")

# Combined plot: shared genes first within each cohort/group, then the remaining genes.
# Sort each block by rho: decreasing for promoting, increasing for suppressive.
immune65_both_first_rho_results <- immune65_results %>%
  mutate(both_studies_first = as.character(gene) %in% c(genes_suppressive_both, genes_promoting_both),
         rho_sort = ifelse(immune_group == "promoting", -rho, rho)) %>%
  arrange(cohort, immune_group, desc(both_studies_first), rho_sort, padj, as.character(gene))
immune65_tcga_both_first_rho_data <- immune65_both_first_rho_results %>%
  filter(cohort == "TCGA") %>%
  mutate(gene = factor(as.character(gene), levels = as.character(gene)))
immune65_cptac_both_first_rho_data <- immune65_both_first_rho_results %>%
  filter(cohort == "CPTAC") %>%
  mutate(gene = factor(as.character(gene), levels = as.character(gene)))

p_tcga_both_first_rho <- p_tcga %+% immune65_tcga_both_first_rho_data
p_cpatc_both_first_rho <- p_cpatc %+% immune65_cptac_both_first_rho_data
p_tcga_cptac_both_first_rho <- patchwork::wrap_plots(
  p_tcga_both_first_rho, p_cpatc_both_first_rho, ncol = 1, guides = "collect"
) & theme(legend.position = "right")
p_tcga_cptac_both_first_rho