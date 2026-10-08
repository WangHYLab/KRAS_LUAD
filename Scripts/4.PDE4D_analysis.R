WORK_DIR<-"path/to/KRAS_LUAD/"
# ============================================================================
# PART 1 — PACKAGES AND HELPER FUNCTIONS
# ============================================================================
# Load required packages
library(DESeq2)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(progeny)
library(ggpubr)
library(dorothea)
library(viper)
preprocess_exp_data<-function(data){
  colnames(data) <- sub(";.*$", "", colnames(data))
  data1<-data[!duplicated(data$gene_name),]
  rownames(data1)<-data1[,'gene_name']
  data2<-data1[,-c(1,2)]
  return(data2)
}

get_deseq2<-function(exp_count,selected_samples){
  countdata <- exp_count[, selected_samples$sample]
  # Build sample metadata
  sample_info <- data.frame(
    row.names = selected_samples$sample,
    group = factor(selected_samples$group, levels = c("low", "high"))
  )
  dds <- DESeqDataSetFromMatrix(countData =countdata,
                                colData = sample_info,
                                design = ~ group)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "high", "low"))
  res <- res[order(res$pvalue), ]
  res_df <- as.data.frame(res)
  res_df$genes<-rownames(res_df)
  return(list(dds=dds,res_df=res_df))
}


# Perform KEGG enrichment analysis
get_KEGG <- function(res, padj_cutoff, log2FC_cutoff, cate) {
  genes_up <- res %>% filter(padj < padj_cutoff, log2FoldChange > log2FC_cutoff) %>% pull(genes)
  genes_down <- res %>% filter(padj < padj_cutoff, log2FoldChange < -log2FC_cutoff) %>% pull(genes)
  
  entrez_up <- if(length(genes_up) > 0) bitr(genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) else data.frame()
  entrez_down <- if(length(genes_down) > 0) bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) else data.frame()
  
  df_up <- data.frame()
  df_down <- data.frame()
  
  # Set a minimum gene-count threshold to prevent enrichKEGG errors
  min_genes <- 5
  
  # KEGG enrichment for upregulated genes
  if (nrow(entrez_up) >= min_genes) {
    kegg_up <- tryCatch({
      enrichKEGG(gene = entrez_up$ENTREZID, organism = "hsa", pAdjustMethod = "BH")
    }, error = function(e) {
      message("KEGG enrichment error (upregulated) in ", cate, ": ", e$message)
      NULL
    })
    if(!is.null(kegg_up)) {
      df_up <- kegg_up@result
      df_up$Group <- "Upregulated"
      df_up$treatment <- cate
      df_up$geneID <- sapply(df_up$geneID, function(ids) {
        genes <- unlist(strsplit(ids, "/"))
        name_df <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        symbol_vec <- name_df$SYMBOL[match(genes, name_df$ENTREZID)]
        paste(symbol_vec, collapse = "/")
      })
      write.csv(df_up, file = paste0('KEGG_results_up_', cate, '.csv'), row.names = FALSE)
    }
  } else {
    message("Too few upregulated genes (", nrow(entrez_up), ") for KEGG in ", cate)
  }
  
  # KEGG enrichment for downregulated genes
  if (nrow(entrez_down) >= min_genes) {
    kegg_down <- tryCatch({
      enrichKEGG(gene = entrez_down$ENTREZID, organism = "hsa", pAdjustMethod = "BH")
    }, error = function(e) {
      message("KEGG enrichment error (downregulated) in ", cate, ": ", e$message)
      NULL
    })
    if(!is.null(kegg_down)) {
      df_down <- kegg_down@result
      df_down$Group <- "Downregulated"
      df_down$treatment <- cate
      df_down$geneID <- sapply(df_down$geneID, function(ids) {
        genes <- unlist(strsplit(ids, "/"))
        name_df <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        symbol_vec <- name_df$SYMBOL[match(genes, name_df$ENTREZID)]
        paste(symbol_vec, collapse = "/")
      })
      write.csv(df_down, file = paste0('KEGG_results_down_', cate, '.csv'), row.names = FALSE)
    }
  } else {
    message("Too few downregulated genes (", nrow(entrez_down), ") for KEGG in ", cate)
  }
  
  return(list(df_up, df_down))
}

# Plot PROGENy pathway activity by group
plot_progeny_box <- function(prog, samples, pathway, group_col = "group") {
  
  df <- data.frame(
    activity = prog[, pathway],
    group = samples[[group_col]]
  )
  y_range <- range(df$activity, na.rm = TRUE)
  y_pos <- y_range[2] + 0.08 * diff(y_range)
  
  ggplot(df, aes(x = group, y = activity, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1) +
    #stat_compare_means(method = "wilcox",comparisons=list(c('low','high')), label = "p.signif") +
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("low", "high")),
      label = "p.format",
      label.y = y_pos
    )+
    theme_classic() +
    labs(
      title = paste("pathway activity:", pathway),
      y = "PROGENy score",
      x=NULL
    )
}

plot_TF_box <- function(tf_scores, samples, tf = "CREB1", group_col = "group") {
  
  df <- data.frame(
    activity = as.numeric(tf_scores[tf, ]),
    group = samples[colnames(tf_scores), group_col]
  )
  y_range <- range(df$activity, na.rm = TRUE)
  y_pos <- y_range[2] + 0.08 * diff(y_range)
  
  ggplot(df, aes(x = group, y = activity, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1) +
    # stat_compare_means(
    #   method = "wilcox",
    #   comparisons = list(c("low", "high")),
    #   label = "p.signif"
    # ) +
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("low", "high")),
      label = "p.format",
      label.y = y_pos
    )+
    theme_classic() +
    labs(
      title = paste("TF activity:", tf),
      y = "VIPER activity score",
      x=NULL
    )
}

# Plot DESeq2-normalized gene expression by group
plot_deseq_gene_boxplot <- function(
    norm_counts,
    coldata,
    res,
    gene_name,
    group_col = "group",
    group_levels=c("high","low"),
    contrast_high = "high",
    contrast_low = "low",
    y_label = "log2(normalized counts + 1)",
    legend_title = "PDE4D"
) {
  
  coldata[[group_col]] <- factor(
    coldata[[group_col]],
    levels = group_levels
  )
  res_gene <- res[gene_name, ]
  
  plot_df <- data.frame(
    sample = rownames(coldata),
    group = coldata[[group_col]],
    expression = log2(norm_counts[gene_name, rownames(coldata)] + 1)
  )
  
  p_label <- paste0("DESeq2 adj. P = ", signif(res_gene$padj, 2))
  
  y_range <- range(plot_df$expression, na.rm = TRUE)
  y_pos <- y_range[2] + 0.08 * diff(y_range)
  
  if (diff(y_range) == 0) {
    y_pos <- y_range[2] + 0.5
  }
  
  p <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 1
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_pos,
      label = p_label,
      size = 3.5
    ) +
    scale_x_discrete(
      limits = group_levels,
      labels = c(
        low = "low",
        high = "high"
      )
    ) +
    labs(
      x = NULL,
      y = y_label,
      title = gene_name,
      fill = legend_title
    ) +
    theme_classic()
  
  return(p)
}

split_by_cutoff <- function(df, expr_col, cutoff, new_col = "group") {
  df[[new_col]] <- ifelse(df[[expr_col]] > cutoff, "high", "low")
  return(df)
}
find_valley_cutoff <- function(x) {
  d <- density(x, na.rm = TRUE)
  y <- d$y
  xgrid <- d$x
  
  # Identify local minima
  dy <- diff(y)
  sign_change <- diff(sign(dy))
  valley_idx <- which(sign_change == 2) + 1
  
  # Select the valley between the two major peaks
  peak_idx <- which(diff(sign(diff(y))) == -2) + 1
  if (length(peak_idx) >= 2 && length(valley_idx) >= 1) {
    # Return the valley located between the first two peaks
    p1 <- peak_idx[1]
    p2 <- peak_idx[2]
    valley_between <- valley_idx[valley_idx > p1 & valley_idx < p2]
    return(xgrid[valley_between][1])
  } else {
    return(NA)
  }
}


# ============================================================================
# PART 2 — PDE4D GROUPING AND DIFFERENTIAL EXPRESSION IN TCGA AND CCLE
# ============================================================================
# TCGA KRAS-mutant LUAD cohort
load(paste0(WORK_DIR,'Rdata/kras_count138.RData'))
load(paste0(WORK_DIR,'Rdata/kras_tpm138.RData'))
count_tcga<-preprocess_exp_data(kras_count138)
tpm_tcga<-preprocess_exp_data(kras_tpm138)
df_tcga <- data.frame(sample = colnames(tpm_tcga),
                 expr = as.numeric(tpm_tcga["PDE4D", ]))
upper_threshold_tcga <- quantile(df_tcga$expr, 0.7)
lower_threshold_tcga <- quantile(df_tcga$expr, 0.3)

df_tcga$group <- NA
df_tcga$group[df_tcga$expr >= upper_threshold_tcga] <- "high"
df_tcga$group[df_tcga$expr <= lower_threshold_tcga] <- "low"

selected_samples_tcga <- df_tcga[!is.na(df_tcga$group), ]

deseq2_tcga<-get_deseq2(exp_count = count_tcga,selected_samples = selected_samples_tcga)
dds_PDE4D_tcga<-deseq2_tcga$dds
save(dds_PDE4D_tcga,file=paste0(WORK_DIR,file='Rdata/dds_PDE4D_tcga.RData'))

res_PDE4D_tcga<-deseq2_tcga$res_df
write.csv(res_PDE4D_tcga,file=paste0(WORK_DIR,file='Results/all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv'),row.names=F)


# CCLE KRAS-mutant LUAD cell lines
ccle_data_dir<-'path/to/CCLE_data/'
# Define PDE4D expression groups
genes_tpm<-read.csv(file=paste0(ccle_data_dir,'DepMap_luad_celllines_summary_with_genes_tpmlogp1.csv'))
kras_luad<-genes_tpm[which(genes_tpm$KRAS_mutation!='wt'),]
# Supplementary Figure 5a: PDE4D expression density
ggplot(kras_luad, aes(x = PDE4D)) +
  geom_density(size = 1) +
  theme_classic() +
  labs(
    title = "Density of PDE4D expression",
    x = "PDE4D expression log2(TPM+1)",
    y = "Density"
  )
cutoff <- find_valley_cutoff(kras_luad$PDE4D)
pde4d_group<-split_by_cutoff(kras_luad,'PDE4D',cutoff = cutoff,new_col = 'group')


# DESeq2 differential expression analysis
rawcount<-read.csv(paste0(ccle_data_dir,'OmicsExpressionRawReadCountHumanAllGenes.csv'))
krasluad_count<-rawcount[which(rawcount$ModelID %in% kras_luad$ModelID),]
colnames(krasluad_count) <- gsub(" ", "", sub("\\.\\.[0-9]+\\.$", "", colnames(krasluad_count)))
krasluad_count<-krasluad_count[which(krasluad_count$IsDefaultEntryForModel=='Yes'),]
rownames(krasluad_count)<-krasluad_count$ModelID

krasluad_count_clean<-krasluad_count[,-c(1:6)]
krasluad_count_clean1<-as.data.frame(t(krasluad_count_clean))

sample_info<-data.frame(sample=pde4d_group$ModelID,
                        group=pde4d_group$group)

sample_info <- sample_info[match(colnames(krasluad_count_clean1), sample_info$sample), ]
rownames(sample_info) <- sample_info$sample
write.csv(sample_info,file = paste0(WORK_DIR,'Results/samples_PDE4D_group_cutoff_valley_ccle.csv'),row.names = F)

# Run DESeq2 and save the results
deseq2_ccle<-get_deseq2(exp_count = krasluad_count_clean1,selected_samples = sample_info)
dds_ccle<-deseq2_ccle$dds
save(dds_ccle,file=paste0(WORK_DIR,'Rdata/dds_PDE4D_ccle.RData'))
res_PDE4D_ccle<-deseq2_ccle$res_df
write.csv(res_PDE4D_ccle,file=paste0(WORK_DIR,'Results/all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv'),row.names = F)


# KEGG enrichment analysis
padj_cutoff<-0.05
log2FC_cutoff<-log2(1.5)
library(openxlsx)

#####CCLE
res_PDE4D_ccle<-read.csv(file=paste0(WORK_DIR,file='Results/all_deseq2_results_KRASluad_PDE4D_high_vs_low_CCLE.csv'))
rownames(res_PDE4D_ccle)<-res_PDE4D_ccle$genes

kegg_res_ccle<-get_KEGG(res = res_PDE4D_ccle,padj_cutoff = padj_cutoff,log2FC_cutoff = log2FC_cutoff,cate = 'PDE4D')
names(kegg_res_ccle)<-c('DEG_up','DEG_down')
lst_safe <- lapply(kegg_res_ccle, function(df) {
  df[] <- lapply(df, as.character)
  df
})
write.xlsx(lst_safe, file=paste0(WORK_DIR,'Results/DEGs_enrichkegg_res_ccle_PDE4D.xlsx'))
save(kegg_res_ccle,file=paste0(WORK_DIR,'Rdata/kegg_res_ccle_PDE4D.RData'))
load(file=paste0(WORK_DIR,'Rdata/kegg_res_ccle_PDE4D.RData'))

#####tcga
res_PDE4D_tcga<-read.csv(file=paste0(WORK_DIR,file='Results/all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv'))
rownames(res_PDE4D_tcga)<-res_PDE4D_tcga$genes

kegg_res_tcga<-get_KEGG(res = res_PDE4D_tcga,padj_cutoff = padj_cutoff,log2FC_cutoff = log2FC_cutoff,cate = 'PDE4D')
names(kegg_res_tcga)<-c('DEG_up','DEG_down')
lst_safe <- lapply(kegg_res_tcga, function(df) {
  df[] <- lapply(df, as.character)
  df
})
library(openxlsx)
write.xlsx(lst_safe, file=paste0(WORK_DIR,'Results/DEGs_enrichkegg_res_tcga_PDE4D.xlsx'))
save(kegg_res_tcga,file=paste0(WORK_DIR,'Rdata/kegg_res_tcga_PDE4D.RData'))

load(file=paste0(WORK_DIR,'Rdata/kegg_res_tcga_PDE4D.RData'))
###kegg_res_tcga <- read.xlsx("Results/DEGs_enrichkegg_res_tcga_PDE4D.xlsx",colTypes = "text")
kegg_res_tcga1 <- lapply(kegg_res_tcga, function(df) {
  rownames(df) <-seq(1:nrow(df)) 
  df
})
kegg_v<-rbind(kegg_res_tcga1$DEG_up[c(1:5),],kegg_res_tcga1$DEG_down[c(1:5),])
p_kegg_tcga<-ggplot(kegg_v, aes(x = GeneRatio, 
                          y = reorder(Description, -log10(pvalue)), 
                          size = Count, 
                          color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    facet_grid(Group ~ ., scales = "free", space = "free") +  # 上下排列
    labs(title = "KRAS-mutant TCGA-LUAD: PDE4D-high vs PDE4D-low",
         x = "Gene Ratio", y = "KEGG pathways", color = "p.adjust") +
    theme_bw() +
    theme(strip.text.y = element_text(size=10, face="bold"),
          axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=5),
          plot.title = element_text(hjust=0.5, size=10),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8))
p_kegg_tcga

# ============================================================================
# PART 3 — IMMUNE-CELL INFILTRATION ANALYSES
# ============================================================================
# ImmuCellAI gene signatures for 24 immune-cell types
sig<-read.csv(paste0(WORK_DIR,'Tables/ImmuCellAI_24_immune_cell_types_signature_genes.csv'))
rownames(sig)<-sig[,1]
sig<-sig[,-1]
sig<-as.matrix(sig)
sigList<-list()# Convert the signature of each cell type to a list
for(i in 1:nrow(sig)){
  x<-as.vector(sig[i,])# Column names can also be removed with names(x) <- NULL
  x<-x[nchar(x)>0]
  sigList[[i]]<-x
}
cell<-rownames(sig)
names(sigList)<-cell# Assign cell-type names to sigList

gene_counts <- sapply(sigList, length)
sig_df <- data.frame(
  term = rep(names(sigList), times = gene_counts),  # Repeat each term according to its gene count
  gene = unlist(sigList, use.names = FALSE)         # Flatten the gene lists without retaining names
)

tcga<-read.csv(paste0(WORK_DIR,'/Results/all_deseq2_PDE4D_high_vs_low_30_tcga_luad_kras.csv'))

stat_tcga <- tcga$stat
names(stat_tcga) <- tcga$genes
stat_tcga <- sort(stat_tcga, decreasing = TRUE)
immune_ges_tcga_stat_pde4d<-GSEA(stat_tcga,TERM2GENE=sig_df,pvalueCutoff = 1)
save(immune_ges_tcga_stat_pde4d,file=paste0(WORK_DIR,'Rdata/immune_ges_tcga_stat_pde4d.Rdata'))
load(file=paste0(WORK_DIR,'Rdata/immune_ges_tcga_stat_pde4d.Rdata'))  
write.csv(immune_ges_tcga_stat_pde4d,file =paste0(WORK_DIR,'Results/immune_ges_tcga_stat_PDE4D.csv'),row.names = F)
immune_res_tcga<-immune_ges_tcga_stat_pde4d@result
immune_res_tcga$annotation<-c(rep('NES<0',20))
immune_res_tcga[which(immune_res_tcga$NES>0),'annotation']<-'NES>0'
colnames(immune_res_tcga)[c(1,4)]<-c('immunecell','ES')
immune_res_tcga$logqvalue<--log10(immune_res_tcga$qvalue)

immune_res_tcga_sorted <- immune_res_tcga %>%
  arrange(logqvalue, desc(Description)) %>%  # Use descending order for Description
  # Convert immunecell to a factor to preserve the sorted plotting order
  mutate(immunecell = factor(immunecell, levels = unique(immunecell)))
# Supplementary Figure 4c: ImmuCellAI GSEA
p_immucellAI<-ggplot(immune_res_tcga_sorted,aes(x=logqvalue,y=immunecell,size=abs(NES),colour=annotation))+
  scale_color_manual(values=c("blue","red"))+
  geom_point(shape=16)+labs(x='-log10(qvalue)',y='Immune cells',title = 'PDE4D high vs PDE4D low (KRAS-mut LUAD)')+
  theme_set(theme_bw())+theme(axis.title = element_text(size=13),axis.title.x = element_text(size=13),legend.title = element_text(size=9))

# CIBERSORT immune-cell analysis
immuneSig<-read.csv(paste0(WORK_DIR,'Tables/CIBERSORT_cell_types_markers.csv'))
# Gene set enrichment analysis
immune_ges_stat_cibersort_pde4d<-GSEA(stat_tcga,TERM2GENE=immuneSig,pvalueCutoff = 1)
save(immune_ges_stat_cibersort_pde4d,file=paste0(WORK_DIR,'Rdata/immune_ges_stat_cibersort_pde4d.RData'))
load(file=paste0(WORK_DIR,'Rdata/immune_ges_stat_cibersort_pde4d.RData'))
write.csv(immune_ges_stat_cibersort_pde4d,file =paste0(WORK_DIR,'Results/immune_ges_stat_cibersort_pde4d.csv'),row.names = F)

# Prepare CIBERSORT GSEA results
immune_res_tcga_cibersort<-immune_ges_stat_cibersort_pde4d@result
immune_res_tcga_cibersort$annotation<-c(rep('NES<0',22))
immune_res_tcga_cibersort[which(immune_res_tcga_cibersort$NES>0),'annotation']<-'NES>0'
colnames(immune_res_tcga_cibersort)[c(1,4)]<-c('immunecell','ES')
immune_res_tcga_cibersort$logqvalue<--log10(immune_res_tcga_cibersort$qvalue)

immune_res_tcga_sorted_cibersort <- immune_res_tcga_cibersort %>%
  arrange(logqvalue, desc(Description)) %>%  # Use descending order for Description
  # Convert immunecell to a factor to preserve the sorted plotting order
  mutate(immunecell = factor(immunecell, levels = unique(immunecell)))
# Supplementary Figure 4b: CIBERSORT GSEA
p_cibersort<-ggplot(immune_res_tcga_sorted_cibersort,aes(x=logqvalue,y=immunecell,size=abs(NES),colour=annotation))+
  scale_color_manual(values=c("blue","red"))+
  geom_point(shape=16)+labs(x='-log10(qvalue)',y='Immune cells cibersort',title = 'PDE4D high vs PDE4D low (KRAS-mut LUAD)')+
  theme_set(theme_bw())+theme(axis.title = element_text(size=13),axis.title.x = element_text(size=13),legend.title = element_text(size=9))



# TIMER2 immune-cell estimation
load(paste0(WORK_DIR,'Rdata/kras_tpm138.RData'))
kras_tpm<-kras_tpm138[!duplicated(kras_tpm138$gene_name),]
kras_tpm_coding<-kras_tpm[which(kras_tpm$genetype=='protein_coding'),]
rownames(kras_tpm_coding)<-kras_tpm_coding$gene_name
kras_exp<-kras_tpm_coding[,-c(1,2)]
write.csv(kras_exp,file=paste0(WORK_DIR,'Results/kras_exp.csv'))

load(file=paste0(WORK_DIR,file='Rdata/dds_PDE4D_tcga.RData'))
samples_pde4d<-as.data.frame(dds_PDE4D_tcga@colData)
kras_exp_pde4d<-kras_exp[,rownames(samples_pde4d)]
# Functions for gene–immune-cell correlation analysis
library(stringr)

get_gene_vs_all_celltypes <- function(exp_mat, cell_mat, gene_name) {
  # Verify that the gene is present
  if (!(gene_name %in% rownames(exp_mat))) {
    stop("gene not found")
  }
  
  # Extract gene expression
  gene_vec <- as.numeric(exp_mat[gene_name, ])
  names(gene_vec) <- colnames(exp_mat)
  
  # Iterate over cell types
  res_list <- lapply(seq_len(nrow(cell_mat)), function(i) {
    cell_name <- cell_mat$cell_type[i]
    cell_vec <- as.numeric(cell_mat[i, -1])
    names(cell_vec) <- colnames(cell_mat)[-1]
    
    # Match samples
    common_samples <- intersect(names(gene_vec), names(cell_vec))
    gene_use <- gene_vec[common_samples]
    cell_use <- cell_vec[common_samples]
    
    # Remove missing values
    keep <- complete.cases(gene_use, cell_use)
    gene_use <- gene_use[keep]
    cell_use <- cell_use[keep]
    
    # Exclude insufficient or invariant data
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
    
    # Calculate the Spearman correlation
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
  
  # Combine results and adjust P values
  res_df <- bind_rows(res_list) %>%
    mutate(
      padj = p.adjust(pvalue, method = "BH")
    ) %>%
    arrange(pvalue)
  
  return(res_df)
}

extract_supported_celltypes <- function(res_df, padj_cutoff = 0.05, min_tools_same_direction = 2) {
  
  res2 <- res_df %>%
    mutate(
      method = str_extract(cell_type, "[^_]+$"),
      cell_type_main = str_remove(cell_type, "_[^_]+$")
    ) %>%
    filter(!is.na(padj), padj < padj_cutoff)
  
  summary_res <- res2 %>%
    group_by(cell_type_main) %>%
    summarise(
      n_sig_tools = n(),
      
      n_positive = sum(direction == "positive", na.rm = TRUE),
      n_negative = sum(direction == "negative", na.rm = TRUE),
      
      positive_tools = paste(method[direction == "positive"], collapse = "; "),
      negative_tools = paste(method[direction == "negative"], collapse = "; "),
      
      positive_celltypes = paste(cell_type[direction == "positive"], collapse = "; "),
      negative_celltypes = paste(cell_type[direction == "negative"], collapse = "; "),
      
      positive_cor = paste(round(cor[direction == "positive"], 3), collapse = "; "),
      negative_cor = paste(round(cor[direction == "negative"], 3), collapse = "; "),
      
      positive_padj = paste(signif(padj[direction == "positive"], 3), collapse = "; "),
      negative_padj = paste(signif(padj[direction == "negative"], 3), collapse = "; "),
      
      mean_cor_positive = ifelse(any(direction == "positive"),
                                 mean(cor[direction == "positive"], na.rm = TRUE), NA),
      mean_cor_negative = ifelse(any(direction == "negative"),
                                 mean(cor[direction == "negative"], na.rm = TRUE), NA),
      .groups = "drop"
    ) %>%
    mutate(
      support_positive = n_positive >= min_tools_same_direction,
      support_negative = n_negative >= min_tools_same_direction,
      supported_direction = case_when(
        support_positive & support_negative ~ "both",
        support_positive ~ "positive",
        support_negative ~ "negative",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(supported_direction)) %>%
    arrange(desc(pmax(n_positive, n_negative)), cell_type_main)
  
  return(summary_res)
}
extract_supported_celltypes_long <- function(res_df, padj_cutoff = 0.05, min_tools_same_direction = 2) {
  
  # Separate cell-type and estimation-method names
  res2 <- res_df %>%
    mutate(
      method = str_extract(cell_type, "[^_]+$"),
      cell_type_main = str_remove(cell_type, "_[^_]+$")
    )
  
  # Count the methods with available data for each cell type
  # A non-missing correlation indicates that the method has available data
  tool_count_df <- res2 %>%
    filter(!is.na(cor)) %>%
    distinct(cell_type_main, method) %>%
    count(cell_type_main, name = "n_available_tools")
  
  # Retain significant results
  sig_res <- res2 %>%
    filter(!is.na(padj), padj < padj_cutoff)
  
  # Summarize positive correlations
  pos_res <- sig_res %>%
    filter(direction == "positive") %>%
    group_by(cell_type_main) %>%
    summarise(
      supported_direction = "positive",
      n_sig_tools = n(),
      tools = paste(method, collapse = "; "),
      full_celltypes = paste(cell_type, collapse = "; "),
      cor_values = paste(round(cor, 3), collapse = "; "),
      pvalue_values = paste(signif(pvalue, 3), collapse = "; "),
      padj_values = paste(signif(padj, 3), collapse = "; "),
      mean_cor = mean(cor, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Summarize negative correlations
  neg_res <- sig_res %>%
    filter(direction == "negative") %>%
    group_by(cell_type_main) %>%
    summarise(
      supported_direction = "negative",
      n_sig_tools = n(),
      tools = paste(method, collapse = "; "),
      full_celltypes = paste(cell_type, collapse = "; "),
      cor_values = paste(round(cor, 3), collapse = "; "),
      pvalue_values = paste(signif(pvalue, 3), collapse = "; "),
      padj_values = paste(signif(padj, 3), collapse = "; "),
      mean_cor = mean(cor, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_res <- bind_rows(pos_res, neg_res) %>%
    left_join(tool_count_df, by = "cell_type_main") %>%
    mutate(
      keep = case_when(
        n_available_tools == 1 & n_sig_tools >= 1 ~ TRUE,
        n_available_tools >= 2 & n_sig_tools >= min_tools_same_direction ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(keep) %>%
    arrange(desc(n_sig_tools), desc(abs(mean_cor)))
  
  return(summary_res)
}

matrix<-read.csv(paste0(WORK_DIR,'Results/TIMER2_estimation_matrix_KRAStpm.csv'))
colnames(matrix) <- gsub("\\.", "-", colnames(matrix))
matrix_pde4d<-matrix[,c("cell_type",rownames(samples_pde4d))]

res_PDE4D <- get_gene_vs_all_celltypes(
  exp_mat = kras_exp,
  cell_mat = matrix,
  gene_name = "PDE4D"
)

# Correlation heatmap
library(tidyr)
res_PDE4D_plot <- res_PDE4D %>%
  mutate(
    cell_type_main = str_remove(cell_type, "_[^_]+$"),
    method = str_extract(cell_type, "[^_]+$"),
    cor = ifelse(is.infinite(cor), NA, cor)
  ) %>%
  group_by(cell_type_main) %>%
  mutate(
    n_tool_valid = n_distinct(method[!is.na(cor)]),
    n_method_p005 = n_distinct(method[pvalue < 0.05 & !is.na(cor)]),
    n_method_p010 = n_distinct(method[pvalue < 0.10 & !is.na(cor)])
  ) %>%
  ungroup()
res_PDE4D_plot_filter <- res_PDE4D_plot %>%
  filter(n_tool_valid > 1)
cell_order <- res_PDE4D_plot_filter %>%
  group_by(cell_type_main) %>%
  summarise(
    n_neg_method = n_distinct(method[cor < 0]),
    mean_cor = mean(cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_neg_method), mean_cor) %>%
  pull(cell_type_main)
plot_df <- res_PDE4D_plot_filter %>%
  dplyr::select(cell_type_main, method, cor, pvalue, padj) %>%
  distinct() %>%
  complete(cell_type_main, method) %>%
  mutate(
    sig_label = case_when(
      !is.na(pvalue) & pvalue < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Supplementary Figure 4a: cross-method correlation heatmap
p <- ggplot(plot_df, aes(x = cell_type_main, y = method, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 4) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey80",
    name = "Spearman\ncorrelation"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p
# Optional alternative summary method
# res_PDE4D_focus <- extract_supported_celltypes(res_PDE4D)
library(tidyverse)

res_PDE4D_focus <- extract_supported_celltypes_long(res_PDE4D)
select_cell_types<-res_PDE4D_focus[c(1:11),]
cell_types_select <- select_cell_types %>%
  pull(full_celltypes) %>%
  str_split(pattern = ";") %>%
  unlist() %>%
  str_trim() %>%
  discard(is_empty) %>%
  unique() %>%
  sort()
res_PDE4D_select<-res_PDE4D[which(res_PDE4D$cell_type %in% cell_types_select),]
plot_df <- res_PDE4D_select %>%
  mutate(
    method = str_extract(cell_type, "[^_]+$"),
    cell_type_main = str_remove(cell_type, "_[^_]+$")
  )

# Count the methods supporting each cell type
tool_count <- plot_df %>%
  group_by(cell_type_main) %>%
  summarise(n_tools = n(), .groups = "drop")

# Merge and sort by method count, followed by correlation
plot_df <- plot_df %>%
  left_join(tool_count, by = "cell_type_main") %>%
  arrange(n_tools, cell_type_main) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

celltype_order <- plot_df %>%
  distinct(cell_type_main, n_tools) %>%
  arrange(desc(n_tools), cell_type_main) %>%
  pull(cell_type_main)

plot_df <- plot_df %>%
  mutate(cell_type_main = factor(cell_type_main, levels = celltype_order))
timer2_spearman_cor_plot_df<-plot_df
save(timer2_spearman_cor_plot_df,file=paste0(WORK_DIR,"/Rdata/timer2_spearman_cor_plot_df.Rdata"))

# Plot gene–cell-type correlations
library(mgcv)
plot_gene_celltype_cor <- function(
    exp_mat,
    cell_mat,
    gene_name,
    cell_name,
    smooth_method = c("gam", "loess", "lm")
) {
  
  smooth_method <- match.arg(smooth_method)
  
  # 1. Extract gene expression
  if (!(gene_name %in% rownames(exp_mat))) {
    stop("gene_name not found in exp_mat")
  }
  gene_vec <- as.numeric(exp_mat[gene_name, ])
  names(gene_vec) <- colnames(exp_mat)
  
  # 2. Extract cell proportion
  if (!(cell_name %in% cell_mat$cell_type)) {
    stop("cell_name not found in cell_mat")
  }
  cell_vec <- as.numeric(cell_mat[cell_mat$cell_type == cell_name, -1])
  names(cell_vec) <- colnames(cell_mat)[-1]
  
  # 3. Match common samples
  common_samples <- intersect(names(gene_vec), names(cell_vec))
  gene_use <- gene_vec[common_samples]
  cell_use <- cell_vec[common_samples]
  
  # 4. Remove NA
  keep <- complete.cases(gene_use, cell_use)
  gene_use <- gene_use[keep]
  cell_use <- cell_use[keep]
  sample_use <- common_samples[keep]
  
  # 5. Spearman correlation
  cor_res <- cor.test(gene_use, cell_use, method = "spearman")
  
  plot_df <- data.frame(
    sample = sample_use,
    expression = gene_use,
    proportion = cell_use
  )
  
  # 6. Select smoother
  smooth_layer <- switch(
    smooth_method,
    "lm" = geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.6,
      color = "black"
    ),
    "loess" = geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.8,
      linewidth = 0.6,
      color = "black"
    ),
    "gam" = geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 4),
      se = FALSE,
      linewidth = 0.6,
      color = "black"
    )
  )
  
  # 7. Plot
  p <- ggplot(plot_df, aes(x = expression, y = proportion)) +
    geom_point(size = 2.5, color = "black", fill = "grey80", shape = 21) +
    smooth_layer +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(gene_name, " vs ", cell_name),
      subtitle = paste0(
        "Spearman rho = ", round(unname(cor_res$estimate), 3),
        ", p = ", signif(cor_res$p.value, 3),
        ", n = ", nrow(plot_df)
      ),
      x = paste0(gene_name, " expression"),
      y = paste0(cell_name, " score")
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold")
    )
  
  return(list(
    plot = p,
    cor_test = cor_res,
    data = plot_df
  ))
}


matrix<-read.csv(paste0(WORK_DIR,'Results/TIMER2_estimation_matrix_KRAStpm.csv'))
kras_exp<-read.csv(file=paste0(WORK_DIR,'Results/kras_exp.csv'),row.names = 1)

# ============================================================================
# PART 4 — FIGURE 5: IMMUNE CORRELATIONS AND CCLE FUNCTIONAL ANALYSES
# ============================================================================
# Figure 5a: PDE4D expression versus immune score
res_plot_immune <- plot_gene_celltype_cor(
  exp_mat = kras_exp,
  cell_mat = matrix,
  gene_name = "PDE4D",
  cell_name = "immune score_XCELL",
  smooth_method = "gam"
)
res_plot_immune$plot

# Figure 5b: selected TIMER2 immune-cell correlations
load(file=paste0(WORK_DIR,"/Rdata/timer2_spearman_cor_plot_df.Rdata"))
plot_df<-timer2_spearman_cor_plot_df[which(timer2_spearman_cor_plot_df$cell_type_main %in% c("T cell CD8+","Myeloid dendritic cell","Macrophage")),]
plot_df$cell_type_main<-factor(plot_df$cell_type_main,levels = c("T cell CD8+", "Myeloid dendritic cell", "Macrophage"))

# Prepare correlation and P-value panels
plot_df2 <- plot_df %>%
  mutate(
    cell_type_main = factor(
      cell_type_main,
      levels = c("T cell CD8+", "Myeloid dendritic cell", "Macrophage")
    ),
    neg_log10_p = -log10(pmax(pvalue, .Machine$double.xmin))
  ) %>%
  arrange(cell_type_main, cor) %>%   # Control the order within each facet
  mutate(
    cell_type = factor(cell_type, levels = rev(unique(cell_type)))
  )
p <- ggplot(plot_df2, aes(x = cor, y = cell_type)) +
  geom_col(
    fill = "grey85",
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_grid(
    cell_type_main ~ .,
    scales = "free_y",
    space = "free"
  ) +
  labs(
    title = "Immune infiltration estimates",
    x = expression("Spearman's " * rho),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_blank(),
    strip.text.y.left = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(t = 8, r = 15, b = 8, l = 15)
  )
p_pval <- ggplot(plot_df2, aes(x = neg_log10_p, y = cell_type)) +
  geom_segment(
    aes(x = 0, xend = neg_log10_p, y = cell_type, yend = cell_type),
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_point(
    size = 2.2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.5
  ) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  facet_grid(
    cell_type_main ~ .,
    scales = "free_y",
    space = "free"
  ) +
  labs(
    x = expression(-log[10](P)),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_blank(),
    strip.text.y.left = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(t = 8, r = 5, b = 8, l = 2)
  )

## CIBERSORT GSEA: show 1/2/3 rows for CD8 T cells, dendritic cells,
## and macrophages, while retaining the same facet heights as p and p_pval.
load(file = paste0(WORK_DIR, "Rdata/immune_ges_stat_cibersort_pde4d.RData"))
cibersort_gsea_res <- immune_ges_stat_cibersort_pde4d@result

cibersort_gsea_key <- tibble::tribble(
  ~ID,                           ~cell_type_main,              ~gsea_cell_type,                 ~gsea_order,
  "T.cells.CD8",                "T cell CD8+",                "T cells CD8",                           1,
  "Dendritic.cells.resting",    "Myeloid dendritic cell",    "Dendritic cells resting",               1,
  "Dendritic.cells.activated",  "Myeloid dendritic cell",    "Dendritic cells activated",             2,
  "Macrophages.M0",             "Macrophage",                "Macrophages M0",                         1,
  "Macrophages.M1",             "Macrophage",                "Macrophages M1",                         2,
  "Macrophages.M2",             "Macrophage",                "Macrophages M2",                         3
)

missing_cibersort_terms <- setdiff(cibersort_gsea_key$ID, cibersort_gsea_res$ID)
if (length(missing_cibersort_terms) > 0) {
  stop(
    "Missing CIBERSORT GSEA terms: ",
    paste(missing_cibersort_terms, collapse = ", ")
  )
}

## Use the existing TIMER2 rows as an invisible y-axis skeleton. This keeps
## all three facet heights exactly parallel to those in p and p_pval.
cibersort_facet_rows <- plot_df2 %>%
  distinct(cell_type_main, cell_type) %>%
  group_by(cell_type_main) %>%
  mutate(facet_row = row_number(), facet_n = n()) %>%
  ungroup()

cibersort_gsea_plot_df <- cibersort_gsea_key %>%
  left_join(
    cibersort_gsea_res %>% dplyr::select(ID, NES, qvalue),
    by = "ID"
  ) %>%
  mutate(
    cell_type_main = factor(
      cell_type_main,
      levels = c("T cell CD8+", "Myeloid dendritic cell", "Macrophage")
    )
  ) %>%
  group_by(cell_type_main) %>%
  arrange(gsea_order, .by_group = TRUE) %>%
  mutate(target_fraction = row_number() / (n() + 1)) %>%
  ungroup() %>%
  left_join(
    cibersort_facet_rows %>% distinct(cell_type_main, facet_n),
    by = "cell_type_main"
  ) %>%
  mutate(facet_row = pmax(1L, pmin(facet_n, round(target_fraction * (facet_n + 1))))) %>%
  left_join(
    cibersort_facet_rows %>% dplyr::select(cell_type_main, facet_row, cell_type),
    by = c("cell_type_main", "facet_row")
  ) %>%
  mutate(gsea_label = gsea_cell_type)

p_gsea <- ggplot() +
  geom_blank(
    data = plot_df2,
    aes(x = 0, y = cell_type)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey50"
  ) +
  geom_point(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main != "Myeloid dendritic cell"),
    aes(x = NES, y = cell_type),
    shape = 23,
    fill = "#6F8FAF",
    color = "black",
    size = 4
  ) +
  geom_point(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main == "Myeloid dendritic cell", gsea_order == 1),
    aes(x = NES, y = cell_type),
    shape = 23,
    fill = "#6F8FAF",
    color = "black",
    size = 4,
    position = position_nudge(y = -0.5)
  ) +
  geom_point(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main == "Myeloid dendritic cell", gsea_order == 2),
    aes(x = NES, y = cell_type),
    shape = 23,
    fill = "#6F8FAF",
    color = "black",
    size = 4,
    position = position_nudge(y = 0.5)
  ) +
  geom_text(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main != "Myeloid dendritic cell"),
    aes(x = NES, y = cell_type, label = gsea_label),
    nudge_x = 0.35,
    hjust = 0,
    size = 2.5,
    lineheight = 0.9
  ) +
  geom_text(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main == "Myeloid dendritic cell", gsea_order == 1),
    aes(x = NES, y = cell_type, label = gsea_label),
    nudge_x = 0.35,
    nudge_y = -0.5,
    hjust = 0,
    size = 2.5,
    lineheight = 0.9
  ) +
  geom_text(
    data = cibersort_gsea_plot_df %>%
      filter(cell_type_main == "Myeloid dendritic cell", gsea_order == 2),
    aes(x = NES, y = cell_type, label = gsea_label),
    nudge_x = 0.35,
    nudge_y = 0.5,
    hjust = 0,
    size = 2.5,
    lineheight = 0.9
  ) +
  coord_cartesian(
    xlim = c(
      min(cibersort_gsea_plot_df$NES, 0, na.rm = TRUE) - 0.3,
      max(cibersort_gsea_plot_df$NES, 0, na.rm = TRUE) + 0.3
    ),
    clip = "off"
  ) +
  facet_grid(
    cell_type_main ~ .,
    scales = "free_y",
    space = "free"
  ) +
  labs(x = "NES", y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 9),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(5.5, 75, 5.5, 2)
  )

p_gsea
# Assemble Figure 5b
p_final_gsea <- p + p_pval + p_gsea +
  plot_layout(widths = c(1.5, 1, 1.25))
p_final_gsea

# Figures 5c–d: CCLE KEGG and cytokine analyses
load(file=paste0(WORK_DIR,'Rdata/kegg_res_ccle_PDE4D.RData'))

# Alternative input: kegg_res_ccle <- read.xlsx("Results/DEGs_enrichkegg_res_ccle_PDE4D.xlsx",colTypes = "text")
kegg_res_ccle1 <- lapply(kegg_res_ccle, function(df) {
  rownames(df) <-seq(1:nrow(df)) 
  df
})
kegg_v_ccle<-rbind(kegg_res_ccle1$DEG_up[c(1:6),],kegg_res_ccle1$DEG_down[c(1:5),])
# Figure 5c: KEGG enrichment dot plot
p_kegg_ccle<-ggplot(kegg_v_ccle, aes(x = GeneRatio, 
                                y = reorder(Description, -log10(pvalue)), 
                                size = Count, 
                                color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(Group ~ ., scales = "free", space = "free") +  # Arrange groups vertically
  labs(title = "DepMap KRAS-mutant LUAD cell lines: PDE4D-high vs PDE4D-low",
       x = "Gene Ratio", y = "KEGG pathways", color = "p.adjust") +
  theme_bw() +
  theme(strip.text.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=5),
        plot.title = element_text(hjust=0.5, size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8))
p_kegg_ccle

ccle_cytokiness_up<-kegg_res_ccle$DEG_up[which(kegg_res_ccle$DEG_up$ID=="hsa04060"),]
ccle_cytokiness_down<-kegg_res_ccle$DEG_down[which(kegg_res_ccle$DEG_down$ID=="hsa04060"),]

samples_ccle<-read.csv(file = paste0(WORK_DIR,'Results/samples_PDE4D_group_cutoff_valley_ccle.csv'))
rownames(samples_ccle)<-samples_ccle$sample
load(file=paste0(WORK_DIR,'/Rdata/dds_PDE4D_ccle.RData'))
norm_counts <- counts(dds_ccle, normalized = TRUE)

plot_cytokines_heatmap_up_down <- function(
    norm_counts,
    samples,
    cytokines_up,
    cytokines_down,
    sample_col = "sample",
    group_col = "group",
    group_levels = c("low", "high"),
    title = "DepMap KRAS-mutant LUAD cell lines"
) {
  up_genes <- unique(unlist(strsplit(cytokines_up$geneID, "/")))
  down_genes <- unique(unlist(strsplit(cytokines_down$geneID, "/")))
  
  gene_info <- dplyr::bind_rows(
    data.frame(gene = up_genes, Direction = "Upregulated"),
    data.frame(gene = down_genes, Direction = "Downregulated")
  ) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(
      Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))
    )
  
  samples_use <- samples
  samples_use[[group_col]] <- factor(samples_use[[group_col]], levels = group_levels)
  samples_use <- samples_use %>%
    dplyr::filter(.data[[sample_col]] %in% colnames(norm_counts)) %>%
    dplyr::arrange(.data[[group_col]], .data[[sample_col]])
  
  genes_use <- intersect(gene_info$gene, rownames(norm_counts))
  missing_genes <- setdiff(gene_info$gene, rownames(norm_counts))
  if (length(missing_genes) > 0) {
    message("Genes not found in norm_counts: ", paste(missing_genes, collapse = ", "))
  }
  
  heat_mat <- log2(norm_counts[genes_use, samples_use[[sample_col]], drop = FALSE] + 1)
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[is.na(heat_mat_z)] <- 0
  
  sample_group_df <- data.frame(
    sample = samples_use[[sample_col]],
    PDE4D_group = factor(
      paste0("PDE4D_", samples_use[[group_col]]),
      levels = c("PDE4D_high", "PDE4D_low")
    )
  )
  
  gene_order <- gene_info %>%
    dplyr::filter(gene %in% genes_use) %>%
    dplyr::arrange(Direction, gene) %>%
    dplyr::pull(gene)
  heat_df <- as.data.frame(heat_mat_z[gene_order, , drop = FALSE]) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "z_score"
    ) %>%
    dplyr::left_join(gene_info, by = "gene") %>%
    dplyr::left_join(sample_group_df, by = "sample")
  
  heat_df$sample <- factor(heat_df$sample, levels = samples_use[[sample_col]])
  heat_df$gene <- factor(heat_df$gene, levels = rev(gene_order))
  heat_df$Direction <- factor(
    heat_df$Direction,
    levels = c("Upregulated", "Downregulated")
  )
  
  ggplot(heat_df, aes(x = sample, y = gene, fill = z_score)) +
    geom_tile() +
    facet_grid(Direction ~ PDE4D_group, scales = "free", space = "free") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Row z-score"
    ) +
    labs(
      title = title,
      x = NULL,
      y = "genes in cytokine-cytokine receptor interaction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
}


load(file = paste0(WORK_DIR,'/Rdata/tpm_ccle1.Rdata'))
up_genes <- unique(unlist(strsplit(ccle_cytokiness_up$geneID, "/")))
down_genes <- unique(unlist(strsplit(ccle_cytokiness_down$geneID, "/")))

## Filter the original cytokine genes by expression prevalence in the CCLE
## samples used in this analysis. tpm_ccle1 is already log2(TPM + 1).
cytokine_samples_ccle <- intersect(samples_ccle$sample, colnames(tpm_ccle1))
if (length(cytokine_samples_ccle) == 0) {
  stop("No samples in samples_ccle were found in tpm_ccle1.")
}
cytokine_genes_ccle <- intersect(
  unique(c(up_genes, down_genes)),
  rownames(tpm_ccle1)
)

cytokine_gene_detection_ccle <- data.frame(
  gene = cytokine_genes_ccle,
  expressed_fraction = rowMeans(
    tpm_ccle1[
      cytokine_genes_ccle,
      cytokine_samples_ccle,
      drop = FALSE
    ] > 0.5,
    na.rm = TRUE
  )
) %>%
  dplyr::mutate(keep = expressed_fraction >= 0.25)

cytokine_genes_ccle_filtered <- cytokine_gene_detection_ccle %>%
  dplyr::filter(keep) %>%
  dplyr::pull(gene)

if (length(cytokine_genes_ccle_filtered) == 0) {
  stop(
    "No cytokine genes passed the filter: log2(TPM + 1) > 0.5 ",
    "in at least 25% of CCLE samples."
  )
}

filter_kegg_gene_ids <- function(kegg_df, genes_keep) {
  kegg_df$geneID <- vapply(
    strsplit(kegg_df$geneID, "/"),
    function(x) paste(intersect(x, genes_keep), collapse = "/"),
    character(1)
  )
  kegg_df
}

ccle_cytokiness_up_filtered <- filter_kegg_gene_ids(
  ccle_cytokiness_up,
  cytokine_genes_ccle_filtered
)
ccle_cytokiness_down_filtered <- filter_kegg_gene_ids(
  ccle_cytokiness_down,
  cytokine_genes_ccle_filtered
)
# Figure 5d: cytokine gene heatmap
p_cytokines_ccle_heatmap1 <- plot_cytokines_heatmap_up_down(
  norm_counts = norm_counts,
  samples = samples_ccle,
  cytokines_up = ccle_cytokiness_up_filtered,
  cytokines_down = ccle_cytokiness_down_filtered
)
p_cytokines_ccle_heatmap1


# CCLE volcano plot inputs
load(file=paste0(WORK_DIR,'/Rdata/dds_PDE4D_ccle.RData'))
load(file = paste0(WORK_DIR,'/Rdata/tpm_ccle1.Rdata'))

# CCLE PDE4D-high versus PDE4D-low volcano plot
volcano_padj_cutoff_ccle <- 0.05
volcano_log2fc_cutoff_ccle <- log2(1.5)
volcano_highlight_genes_ccle <- c(
  "IL23A", "CSF2", "GDF15", "CCL2", "TGFB3", "CXCL8"
)

volcano_res_ccle <- DESeq2::results(
  dds_ccle,
  contrast = c("group", "high", "low")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      !is.na(padj) & padj < volcano_padj_cutoff_ccle &
        log2FoldChange > volcano_log2fc_cutoff_ccle ~ "Upregulated",
      !is.na(padj) & padj < volcano_padj_cutoff_ccle &
        log2FoldChange < -volcano_log2fc_cutoff_ccle ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    significance = factor(
      significance,
      levels = c("Downregulated", "Not significant", "Upregulated")
    ),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    highlight = gene %in% volcano_highlight_genes_ccle
  )

missing_volcano_highlights_ccle <- setdiff(
  volcano_highlight_genes_ccle,
  volcano_res_ccle$gene
)
if (length(missing_volcano_highlights_ccle) > 0) {
  message(
    "Highlighted genes not found in DESeq2 results: ",
    paste(missing_volcano_highlights_ccle, collapse = ", ")
  )
}

volcano_highlight_df_ccle <- volcano_res_ccle %>%
  dplyr::filter(highlight, is.finite(log2FoldChange), is.finite(neg_log10_padj))
# Supplementary Figure 5b: CCLE volcano plot
p_volcano_ccle <- ggplot(
  volcano_res_ccle %>%
    dplyr::filter(is.finite(log2FoldChange), is.finite(neg_log10_padj)),
  aes(x = log2FoldChange, y = neg_log10_padj)
) +
  geom_point(
    aes(color = significance),
    size = 1.4,
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = c(-volcano_log2fc_cutoff_ccle, volcano_log2fc_cutoff_ccle),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(volcano_padj_cutoff_ccle),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    data = volcano_highlight_df_ccle,
    shape = 21,
    size = 2.8,
    stroke = 0.6,
    fill = "#FFD54F",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "#2166AC",
      "Not significant" = "grey75",
      "Upregulated" = "#B2182B"
    ),
    drop = FALSE
  ) +
  labs(
    title = "DepMap KRAS-mutant LUAD cell lines",
    subtitle="PDE4D-high vs PDE4D-low",
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " adjusted P"),
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano_ccle <- p_volcano_ccle +
    ggrepel::geom_text_repel(
      data = volcano_highlight_df_ccle,
      aes(label = gene),
      size = 3,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 123
    )
} else {
  p_volcano_ccle <- p_volcano_ccle +
    geom_text(
      data = volcano_highlight_df_ccle,
      aes(label = gene),
      nudge_y = 0.25,
      size = 3,
      check_overlap = TRUE
    )
}

p_volcano_ccle

# ============================================================================
# PART 5 — FIGURE 6: ANTIGEN PRESENTATION AND SIGNALING ACTIVITY
# ============================================================================
# Figure 6e: DESeq2-normalized CCLE expression
load(file=paste0(WORK_DIR,'/Rdata/dds_PDE4D_ccle.RData'))
norm_counts <- counts(dds_ccle, normalized = TRUE)
coldata<-dds_ccle@colData
res <- results(
  dds_ccle,
  contrast = c("group", "high", "low")
)
# Figure 6e: HLA gene expression
p_hla_a<-plot_deseq_gene_boxplot(norm_counts = norm_counts,coldata = coldata,res = res,gene_name = 'HLA.A')
p_hla_b<-plot_deseq_gene_boxplot(norm_counts = norm_counts,coldata = coldata,res = res,gene_name = 'HLA.B')
p_hla_c<-plot_deseq_gene_boxplot(norm_counts = norm_counts,coldata = coldata,res = res,gene_name = 'HLA.C')
# Figure 6f: CREB1 and NLRC5 expression
p_creb1<-plot_deseq_gene_boxplot(norm_counts = norm_counts,coldata = coldata,res = res,gene_name = 'CREB1')

p_nlrc5<-plot_deseq_gene_boxplot(norm_counts = norm_counts,coldata = coldata,res = res,gene_name = 'NLRC5')


# Infer CCLE pathway activity with PROGENy
samples_ccle<-read.csv(file = paste0(WORK_DIR,'Results/samples_PDE4D_group_cutoff_valley_ccle.csv'))
rownames(samples_ccle)<-samples_ccle$sample
tpm<-read.csv('xxDepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv')
colnames(tpm)[1]<-'ModelID'
colnames(tpm) <- gsub(" ", "", sub("\\.\\.[0-9]+\\.$", "", colnames(tpm)))

tpm_ccle <- tpm[match(samples_ccle$sample, tpm$ModelID), ]
rownames(tpm_ccle)<-tpm_ccle$ModelID
tpm_ccle1<-t(tpm_ccle[,-1])

load(file=paste0(WORK_DIR,"/Rdata/tpm_ccle1.Rdata"))

prog_ccle <- progeny(
  tpm_ccle1,
  scale = TRUE,
  organism = "Human",
  top = 100
)
all(rownames(prog_ccle) == rownames(samples_ccle)) # Expected: TRUE

p_JAK_STAT_ccle<-plot_progeny_box(prog_ccle, samples_ccle, "JAK-STAT")

# Infer CCLE transcription-factor activity with DoRothEA and VIPER
colnames(dorothea_hs)
regulon <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
regulon_list <- dorothea::df2regulon(regulon)

viper_scores_ccle <- viper(tpm_ccle1, regulon_list, method = "scale")
save(viper_scores_ccle,file=paste0(WORK_DIR,'Rdata/viper_scores_PDE4D_ccle.RData'))

load(file=paste0(WORK_DIR,'Rdata/viper_scores_PDE4D_ccle.RData'))

# JAK–STAT-related activity
p_STAT1<-plot_TF_box(viper_scores_ccle, samples_ccle, tf = "STAT1")

library(gridExtra)
# Figure 6g: combined CCLE panels
grid.arrange(p_JAK_STAT_ccle,p_STAT1,p_nlrc5,ncol=3)


# Supplementary figure: TCGA expression panels
load(file=paste0(WORK_DIR,'/Rdata/dds_PDE4D_tcga.RData'))
norm_counts_tcga <- counts(dds_PDE4D_tcga, normalized = TRUE)
coldata_tcga<-dds_PDE4D_tcga@colData
res_tcga <- results(
  dds_PDE4D_tcga,
  contrast = c("group", "high", "low")
)

samples_tcga<-as.data.frame(dds_PDE4D_tcga@colData)
samples_tcga$group <- factor(samples_tcga$group, levels = c("high", "low"))

p_hla_a_tcga<-plot_deseq_gene_boxplot(norm_counts = norm_counts_tcga,coldata = coldata_tcga,res = res_tcga,gene_name = 'HLA-A')
p_hla_b_tcga<-plot_deseq_gene_boxplot(norm_counts = norm_counts_tcga,coldata = coldata_tcga,res = res_tcga,gene_name = 'HLA-B')
p_hla_c_tcga<-plot_deseq_gene_boxplot(norm_counts = norm_counts_tcga,coldata = coldata_tcga,res = res_tcga,gene_name = 'HLA-C')
p_creb1_tcga<-plot_deseq_gene_boxplot(norm_counts = norm_counts_tcga,coldata = coldata_tcga,res = res_tcga,gene_name = 'CREB1')
p_nlrc5_tcga<-plot_deseq_gene_boxplot(norm_counts = norm_counts_tcga,coldata = coldata_tcga,res = res_tcga,gene_name = 'NLRC5')

# Infer TCGA pathway activity with PROGENy
load(file=paste0(WORK_DIR,file='Rdata/dds_PDE4D_tcga.RData'))
load(paste0(WORK_DIR,'Rdata/kras_tpm138.RData'))
tpm_tcga<-as.matrix(preprocess_exp_data(kras_tpm138))
samples_tcga<-as.data.frame(dds_PDE4D_tcga@colData)
samples_tcga$group <- factor(samples_tcga$group, levels = c("high", "low"))


tpm_tcga1<-as.matrix(tpm_tcga[,which(colnames(tpm_tcga) %in% rownames(samples_tcga))])
prog_tcga <- progeny(
  tpm_tcga1,
  scale = TRUE,
  organism = "Human",
  top = 100
)
save(prog_tcga,file=paste0(WORK_DIR,'Rdata/prog_PDE4D_30_tcga.RData'))

load(file=paste0(WORK_DIR,'Rdata/prog_PDE4D_30_tcga.RData'))
p_JAK_STAT_tcga<-plot_progeny_box(prog_tcga, samples_tcga, "JAK-STAT")

# Infer TCGA transcription-factor activity with DoRothEA and VIPER
colnames(dorothea_hs)
regulon <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
regulon_list <- dorothea::df2regulon(regulon)

viper_scores_tcga <- viper(tpm_tcga1, regulon_list, method = "scale")
save(viper_scores_tcga,file=paste0(WORK_DIR,'Rdata/viper_scores_PDE4D_30_tcga.RData'))

load(file=paste0(WORK_DIR,'Rdata/viper_scores_PDE4D_30_tcga.RData'))
p_STAT1_tf_tcga<-plot_TF_box(viper_scores_tcga, samples_tcga, tf = "STAT1")

grid.arrange(p_hla_a_tcga,p_hla_b_tcga,p_hla_c_tcga,p_JAK_STAT_tcga,p_STAT1_tf_tcga,p_creb1_tcga,ncol=3)

# ============================================================================
# PART 6 — SUPPLEMENTARY FIGURE 3: cAMP-RELATED EXPRESSION HEATMAPS
# ============================================================================
library(tibble)
library(viridis)
load(paste0(WORK_DIR,'Rdata/kras_tpm138.RData'))
# CIIR: cancer-cell-intrinsic immune-related genes in KRAS-mutant LUAD.
ciir<-read.csv(paste0(WORK_DIR,'/Tables/immune_genes_basical_info.csv'))
ciir1<-ciir[,c(1,14)]
ciir_tpm<-kras_tpm138[which(kras_tpm138$gene_name %in% ciir$gene_name),]
rownames(ciir_tpm)<-ciir_tpm$gene_name
ciir_tpm1<-ciir_tpm[,-c(1,2)]
log_tpm <- log2(ciir_tpm1 + 1)

camp_genes<-c("HCAR2","PDE4D","HCAR3","GRIN1","CGA","NPY")
camp_log_tpm<-log_tpm[which(rownames(log_tpm) %in% camp_genes),]

# Order genes by mean log2(TPM + 1)
gene_order <- rowMeans(camp_log_tpm, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

plot_long <- camp_log_tpm[gene_order, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_order)  )
# Supplementary Figure 3a: TCGA heatmap
p_camp_tcga<-ggplot(plot_long, aes(x = gene, y = sample, fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "log2(TPM + 1)",
    option = "viridis"
  ) +
  labs(
    x = NULL,
    y = "KRAS-mutant TCGA-LUAD tumors",
    title = NULL
    # subtitle = paste0(
    #   length(gene_order),
    #   " genes ordered by mean log2(TPM + 1)"
    # )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Prepare CCLE expression data
samples_ccle<-read.csv(file = paste0(WORK_DIR,'Results/samples_PDE4D_group_cutoff_valley_ccle.csv'))
rownames(samples_ccle)<-samples_ccle$sample

tpm<-read.csv('path/to/DepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv')
colnames(tpm)[1]<-'ModelID'
colnames(tpm) <- gsub(" ", "", sub("\\.\\.[0-9]+\\.$", "", colnames(tpm)))

tpm_ccle <- tpm[match(samples_ccle$sample, tpm$ModelID), ]
rownames(tpm_ccle)<-tpm_ccle$ModelID
tpm_ccle1<-t(tpm_ccle[,-1])
load(file = paste0(WORK_DIR,'/Rdata/tpm_ccle1.Rdata'))


camp_tpm_ccle<-tpm_ccle1[which(rownames(tpm_ccle1) %in% camp_genes),]
# Order genes by mean log2(TPM + 1)
gene_order_ccle <- rowMeans(camp_tpm_ccle, na.rm = TRUE) %>%
  sort(decreasing = TRUE) %>%
  names()

plot_long_ccle <- camp_tpm_ccle[gene_order_ccle, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_order_ccle)  )
# Supplementary Figure 3b: CCLE heatmap
p_camp_ccle<-ggplot(plot_long_ccle, aes(x = gene, y = sample, fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "log2(TPM + 1)",
    option = "viridis"
  ) +
  labs(
    x = NULL,
    y = "DepMap KRAS-mutant LUAD cell lines",
    title = NULL
    # subtitle = paste0(
    #   length(gene_order),
    #   " genes ordered by mean log2(TPM + 1)"
    # )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )
grid.arrange(p_camp_tcga,p_camp_ccle,ncol=2)

