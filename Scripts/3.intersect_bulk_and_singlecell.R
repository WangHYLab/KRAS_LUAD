# Figure overview:
#   Fig. 4b: Bulk RNA-seq expression heatmap
#   Fig. 4c: Single-cell expression heatmap
#   Fig. 4d: KEGG enrichment analysis
#   Fig. 4e: Immune-suppressive CIIR genes after 24-h KRAS inhibition
#   Fig. 4f: Immune-suppressive CIIR genes after 72-h KRAS inhibition
#   Fig. 4g: Immune-promoting CIIR genes after 24-h KRAS inhibition
#   Fig. 4h: Immune-promoting CIIR genes after 72-h KRAS inhibition

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
sc<-read.csv('sc_result_pearson_immun_malignantex_normalize_scale_allgenes.csvresult_pearson_immun_malignantex_normalize_scale_nonecopykat_0626_allgenes.csv')
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
write.csv(res_df_72h1,file='./Results/GSE275534_adag_72h_vs_control_deseq2.csv',row.names = F)


# Prepare GSE275534 differential-expression gene sets.
GSE275534_h358_72h<-read.csv(paste0(WORK_DIR,"Results/GSE275534_h358_adag_72h_vs_control_deseq2.csv"))
GSE275534_h358_72h<-GSE275534_h358_72h[which(GSE275534_h358_72h$hgnc_symbol!=''),]
GSE275534_h358_72h_deg<-get_deg(GSE275534_h358_72h,gene_col = 'hgnc_symbol',lfc_cut=0)
# Load immune-promoting and immune-suppressive CIIR genes.
ciir_genes<-read.csv(paste0(WORK_DIR,"Results/fig4.overlap_65_genes_of_bulk_and_sc.csv"))
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
########Fig4g#####
p1 <- plot_log2FC_single(
  h358_24h_genes,
  group_to_plot = "promoting",
  group_label = "immune-promoting CIIR genes",
  fill_color = "grey80",
  border_color = "black"
)
########Fig4e#####
p2 <- plot_log2FC_single(
  h358_24h_genes,
  group_to_plot = "suppressive",
  group_label = "immune-suppressive CIIR genes",
  #fill_color = "grey35"
  fill_color = "grey80",
  border_color = "black"
)

########Fig4h#####
p3 <- plot_log2FC_single(
  h358_72h_genes,
  group_to_plot = "promoting",
  group_label = "immune-promoting CIIR genes",
  fill_color = "grey80",
  border_color = "black"
)
########Fig4f#####
p4 <- plot_log2FC_single(
  h358_72h_genes,
  group_to_plot = "suppressive",
  group_label = "immune-suppressive CIIR genes",
  #fill_color = "grey35"
  fill_color = "grey80",
  border_color = "black"
)
grid.arrange(p2,p4,p1,p3,ncol=4,widths=c(4,2,2,1.5))
