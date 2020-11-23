### Correlation between methylation and gene expression

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(ggpubr)

### projects with > 150 age-DEGs and > 150 age-DMGs
projects <- c("LGG", "BRCA", "UCEC", "ESCA", "KIRP", "OV", "LIHC", "LAML", "SKCM", "PRAD")

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to test correlation for each gene
corr_methy_expr_gene <- function(gene, df_methy, df_expr){
  print(gene)
  methy <- t(df_methy[rownames(df_methy) == gene,])
  colnames(methy) <- "methylation"
  expr <- t(df_expr[rownames(df_expr) == gene,])
  colnames(expr) <- "expression"
  tmp_df <- merge(methy, expr, by = "row.names")
  cor_result <- cor.test(tmp_df$methylation, tmp_df$expression, method = "pearson")
  names(cor_result$estimate) <- gene
  res <- as.data.frame(cor_result$estimate)
  colnames(res) <- "Pearson_cor"
  return(res)
}

### function to quantify correlation for a cancer
corr_methy_expr <- function(project){
  
  print(paste0("Working on: ", project))
  # read files
  methylation <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  DMGs <- as.character(methylation[methylation$Sig == TRUE,]$gene)
  
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  DEGs <- as.character(expression[expression$Sig == TRUE,]$gene)

  overlap_genes <- DMGs[DMGs %in% DEGs]
  
  DMG_only <- DMGs[!(DMGs %in% DEGs)]
  DEG_only <- DEGs[!(DEGs %in% DMGs)]
  
  all_genes <- as.character(expression$gene)
  other_genes <- all_genes[!(all_genes %in% c(overlap_genes, DMG_only, DEG_only))]
  
  ### read methylation data
  df_methy <- read.csv(paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", collapse = ""))
  df_methy <- df_methy[df_methy$X %in% all_genes,]  # only protein-coding genes in biomaRt
  rownames(df_methy) <- df_methy$X
  df_methy$X <- NULL
  colnames(df_methy) <- unlist(lapply(colnames(df_methy), clean_id))
  
  ### read expression data
  df_expr <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered.csv", collapse = ""))
  df_expr <- df_expr[df_expr$gene %in% all_genes,]
  df_expr <- aggregate(df_expr[, -c(1)],
                       by = list(gene = df_expr$gene),
                       FUN = mean,
                       na.rm = TRUE)
  rownames(df_expr) <- df_expr$gene
  df_expr$gene <- NULL
  colnames(df_expr) <- unlist(lapply(colnames(df_expr), clean_id))
  dim(df_expr)
  
  df_methy <- df_methy[,colnames(df_methy) %in% colnames(df_expr)]  # select only samples in both methylation and gene expression
  df_expr <- df_expr[,colnames(df_expr) %in% colnames(df_methy)]    # select only samples in both methylation and gene expression
  df_expr <- log2(df_expr + 1)
  
  # test correlation for each set of genes
  result_overlap_genes <- lapply(overlap_genes, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_DMG_only <- lapply(DMG_only, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_DEG_only <- lapply(DEG_only, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_other <- lapply(other_genes, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  
  result_overlap_genes <- do.call(rbind, result_overlap_genes)
  result_DMG_only <- do.call(rbind, result_DMG_only)
  result_DEG_only <- do.call(rbind, result_DEG_only)
  result_other <- do.call(rbind, result_other)
  
  df_overlap <- as.data.frame(cbind(condition = rep("age_DMGs_DEGs", nrow(result_overlap_genes)), result_overlap_genes))
  df_DMG_only <- as.data.frame(cbind(condition = rep("age_DMGs", nrow(result_DMG_only)), result_DMG_only))
  df_DEG_only <- as.data.frame(cbind(condition = rep("age_DEGs", nrow(result_DEG_only)), result_DEG_only))
  df_other <- as.data.frame(cbind(condition = rep("others", nrow(result_other)), result_other))
  
  df_all <- rbind(df_overlap, df_DMG_only, df_DEG_only, df_other)
  
  #df_all$corr <- as.numeric(as.character(df_all$Pearson_cor))
  write.csv(df_all, paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_corr_coeff_methy_expr.csv", collapse = ""), row.names = TRUE)
  
  ### Supplementary Fig. 14b
  # write source data
  write.csv(df_all, paste0("Source_Data/Supplementary_Fig_14b_", project, ".csv", collapse = ""), row.names = TRUE)
  
  ### density plot
  pdf(paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_corr_coeff_methy_expr.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggplot(df_all, aes(x=Pearson_cor, fill=condition)) + geom_density(alpha = 0.3) +
    scale_fill_manual(values=c("#bd0026", "#225ea8", "#238b45", "#feb24c")) +
    ggtitle(project) +
    xlab("Pearson correlation coefficient") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_text(size = 12,face="bold"),
          legend.text = element_text(size = 12),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  print(p)
  dev.off()
}

lapply(projects[projects!="UCEC"], corr_methy_expr)   # because UCEC has two RNA-seq platforms


### function to quantify correlation for UCEC
corr_methy_expr_UCEC <- function(project){
  
  print(paste0("Working on: ", project))
  # read files
  methylation <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  DMGs <- as.character(methylation[methylation$Sig == TRUE,]$gene)
  
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  DEGs <- as.character(expression[expression$Sig == TRUE,]$gene)
  
  overlap_genes <- DMGs[DMGs %in% DEGs]
  
  DMG_only <- DMGs[!(DMGs %in% DEGs)]
  DEG_only <- DEGs[!(DEGs %in% DMGs)]
  
  all_genes <- as.character(expression$gene)
  other_genes <- all_genes[!(all_genes %in% c(overlap_genes, DMG_only, DEG_only))]
  
  ### read methylation data
  df_methy <- read.csv(paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", collapse = ""))
  df_methy <- df_methy[df_methy$X %in% all_genes,]  # only protein-coding genes in biomaRt
  rownames(df_methy) <- df_methy$X
  df_methy$X <- NULL
  colnames(df_methy) <- unlist(lapply(colnames(df_methy), clean_id))
  
  ### read expression data
  df_expr <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered.csv", collapse = ""))
  df_expr <- df_expr[df_expr$gene %in% all_genes,]
  df_expr <- aggregate(df_expr[, -c(1)],
                       by = list(gene = df_expr$gene),
                       FUN = mean,
                       na.rm = TRUE)
  rownames(df_expr) <- df_expr$gene
  df_expr$gene <- NULL
  colnames(df_expr) <- unlist(lapply(colnames(df_expr), clean_id))
  dim(df_expr)
  
  df_expr_GA <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered_GA_not_overlap_with_Hiseq.csv", collapse = ""))
  df_expr_GA <- df_expr_GA[df_expr_GA$gene %in% all_genes,]
  df_expr_GA <- aggregate(df_expr_GA[, -c(1)],
                       by = list(gene = df_expr_GA$gene),
                       FUN = mean,
                       na.rm = TRUE)
  rownames(df_expr_GA) <- df_expr_GA$gene
  df_expr_GA$gene <- NULL
  colnames(df_expr_GA) <- unlist(lapply(colnames(df_expr_GA), clean_id))
  dim(df_expr_GA)
  
  # join expression from 2 platforms
  df_expr <- cbind(df_expr, df_expr_GA)
  
  df_methy <- df_methy[,colnames(df_methy) %in% colnames(df_expr)]  # select only samples in both methylation and gene expression
  df_expr <- df_expr[,colnames(df_expr) %in% colnames(df_methy)]    # select only samples in both methylation and gene expression
  df_expr <- log2(df_expr + 1)
  
  # test correlation for each set of genes
  result_overlap_genes <- lapply(overlap_genes, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_DMG_only <- lapply(DMG_only, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_DEG_only <- lapply(DEG_only, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  result_other <- lapply(other_genes, corr_methy_expr_gene, df_methy = df_methy, df_expr = df_expr)
  
  result_overlap_genes <- do.call(rbind, result_overlap_genes)
  result_DMG_only <- do.call(rbind, result_DMG_only)
  result_DEG_only <- do.call(rbind, result_DEG_only)
  result_other <- do.call(rbind, result_other)
  
  df_overlap <- as.data.frame(cbind(condition = rep("age_DMGs_DEGs", nrow(result_overlap_genes)), result_overlap_genes))
  df_DMG_only <- as.data.frame(cbind(condition = rep("age_DMGs", nrow(result_DMG_only)), result_DMG_only))
  df_DEG_only <- as.data.frame(cbind(condition = rep("age_DEGs", nrow(result_DEG_only)), result_DEG_only))
  df_other <- as.data.frame(cbind(condition = rep("others", nrow(result_other)), result_other))
  
  df_all <- rbind(df_overlap, df_DMG_only, df_DEG_only, df_other)
  #df_all$corr <- as.numeric(as.character(df_all$corr))
  write.csv(df_all, paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_corr_coeff_methy_expr.csv", collapse = ""), row.names = TRUE)
  
  # write source data
  write.csv(df_all, paste0("Source_Data/Supplementary_Fig_14b_", project, ".csv", collapse = ""), row.names = TRUE)
  
  ### density plot
  pdf(paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_corr_coeff_methy_expr.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggplot(df_all, aes(x=Pearson_cor, fill=condition)) + geom_density(alpha = 0.3) +
    scale_fill_manual(values=c("#bd0026", "#225ea8", "#238b45", "#feb24c")) +
    ggtitle(project) +
    xlab("Pearson correlation coefficient") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_text(size = 14,face="bold"),
          legend.text = element_text(size = 14),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  print(p)
  dev.off()

}

corr_methy_expr_UCEC("UCEC")
