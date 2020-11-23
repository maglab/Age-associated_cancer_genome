### GSEA enrichment using clusterProfiler

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(ggplot2)

### projects with > 150 age-DEGs and > 150 age-DMGs
projects <- c("LGG", "BRCA", "UCEC", "ESCA", "KIRP", "OV", "LIHC", "LAML", "SKCM", "PRAD")

# biomaRt
ensembl_100 <- useMart(host='http://apr2020.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
bmIDs <- getBM(attributes=c('external_gene_name', 'entrezgene_id'), mart = ensembl_100)

### Function for GSEA
GSEA_GO <- function(project, bmIDs){
  
  print(paste0("Start: ", project))
  
  ### Expression ##################################################################################
  df_expr <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", 
                             project, "_gene_expression_with_age.csv", collapse = ""))
  df_expr <- merge(df_expr, bmIDs, by.x = "gene", by.y = "external_gene_name")
  df_expr <- df_expr[complete.cases(df_expr),]
  
  geneList_expr <- df_expr[,"statistic"]
  names(geneList_expr) <- as.character(df_expr[,"entrezgene_id"])
  geneList_expr <- sort(geneList_expr, decreasing = TRUE)
  
  ### GSEA GO
  gsea_GO_expr <- gseGO(geneList     = geneList_expr,
                        OrgDb        = org.Hs.eg.db,
                        ont          = "all",
                        nPerm        = 1000,
                        minGSSize    = 100,
                        maxGSSize    = 500,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
  result_GO_expr <- gsea_GO_expr@result
  write.csv(result_GO_expr, paste0("Analysis_results/Methylation/5_GSEA/", project, "_GSEA_GO_expression.csv", collapse = ""), row.names = FALSE)
  
  ### Methylation #################################################################################
  df_methy <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", 
                              project, "_methylation_with_age.csv", collapse = ""))
  df_methy <- merge(df_methy, bmIDs, by.x = "gene", by.y = "external_gene_name")
  df_methy <- df_methy[complete.cases(df_methy),]
  
  geneList_methy<- df_methy[,"statistic"]
  names(geneList_methy) <- as.character(df_methy[,"entrezgene_id"])
  geneList_methy <- sort(geneList_methy, decreasing = TRUE)
  
  ### GSEA GO
  gsea_GO_methy <- gseGO(geneList     = geneList_methy,
                         OrgDb        = org.Hs.eg.db,
                         ont          = "all",
                         nPerm        = 1000,
                         minGSSize    = 100,
                         maxGSSize    = 500,
                         pvalueCutoff = 0.1,
                         verbose      = FALSE)
  result_GO_methy <- gsea_GO_methy@result
  write.csv(result_GO_methy, paste0("Analysis_results/Methylation/5_GSEA/", project, "_GSEA_GO_methylation.csv", collapse = ""), row.names = FALSE)
  
  result_GO_expr$group <- rep("Expression", nrow(result_GO_expr))
  result_GO_methy$group <- rep("Methylation", nrow(result_GO_methy))
  
  plot_df <- rbind(result_GO_methy,result_GO_expr)
  shared_ID <- names(table(plot_df$ID)[table(plot_df$ID) > 1])  # GO ID that presented in both methylation and expression
  the_rest <- plot_df$ID[!(plot_df$ID %in% shared_ID)]
  
  shared_df <- plot_df[plot_df$ID %in% shared_ID,]
  shared_df <- shared_df[order(shared_df$pvalue),]
  
  the_rest_df <- plot_df[plot_df$ID %in% the_rest,]
  the_rest_df <- the_rest_df[order(the_rest_df$pvalue),]
  
  merged_df <- rbind(shared_df, the_rest_df)
  merged_df$direction <- ifelse(merged_df$enrichmentScore > 0, "Up", "Down")
  merged_df$group_1 <- paste(merged_df$group, merged_df$direction)
  
  table(merged_df$group_1)
  
  # keep only top 8 for each group for plotting
  ID_1 <- as.character(merged_df[merged_df$group_1 == "Expression Down",]$ID[1:8])
  ID_2 <- as.character(merged_df[merged_df$group_1 == "Expression Up",]$ID[1:8])
  ID_3 <- as.character(merged_df[merged_df$group_1 == "Methylation Down",]$ID[1:8])
  ID_4 <- as.character(merged_df[merged_df$group_1 == "Methylation Up",]$ID[1:8])
  
  all_IDs <- c(ID_1, ID_2, ID_3, ID_4)
  all_IDs <- all_IDs[!(is.na(all_IDs))]
  merged_df_1 <- merged_df[merged_df$ID %in% all_IDs,]
  merged_df_1 <- merged_df_1[unlist(lapply(as.character(merged_df_1$Description), nchar)) < 60,]  # remove too long GO term
  
  ### Plot
  ### Fig. 6e and Supplementary Fig. 11
  pdf(paste0("Analysis_results/Methylation/5_GSEA/", project, "_GSEA_GO_expr_methy.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE) 
  p <- ggplot(merged_df_1, aes(x = group, y = Description)) + 
    geom_point(aes(fill = enrichmentScore, size = -log(p.adjust)), colour="black", pch=21) +
    scale_size_continuous(range = c(1, 4)) +
    scale_fill_gradient2(mid = "#bdbdbd", midpoint = 0, low = "#1D91C0", high = "#a50f15", guide = "colourbar") +
    ggtitle(paste0("GSEA: ", project, collapse = "")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 9.5),
          axis.text.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          panel.border = element_rect(linetype = "solid", fill = NA),
          panel.background = element_blank(),
          #axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "#d9d9d9"),
          panel.grid.minor = element_line(colour = "#d9d9d9"))
  print(p)
  dev.off()

  print(paste0("Finish: ", project))
}

lapply(projects, GSEA_GO, bmIDs = bmIDs)

