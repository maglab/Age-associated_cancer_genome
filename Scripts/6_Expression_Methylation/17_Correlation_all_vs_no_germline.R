### Correlation between fold changes with age of gene expression from all samples and from samples without germline variants

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)

test_corr <- function(project){
  
  # read result from all samples
  All <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  All_coeff <- All[,c("gene", "estimate")]
  colnames(All_coeff) <- c("gene", "All_samples")
  
  # read result from samples without germline mutations
  no_germline <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age_exclude_germline.csv", collapse = ""))
  no_germline_coeff <- no_germline[,c("gene", "estimate")]
  colnames(no_germline_coeff) <- c("gene", "No_germline")
  
  # merge
  df <- merge(All_coeff, no_germline_coeff, by.x = "gene", by.y = "gene")
  
  # write source data
  write.csv(df, paste0("Source_Data/Supplementary_Fig_12a_expression_", project, ".csv", collapse = ""), row.names = FALSE)
  
  # plot
  xmax <- max(abs(min(df$All_samples)), abs(max(df$All_samples)))
  ymax <- max(abs(min(df$No_germline)), abs(max(df$No_germline)))
  
  pdf(paste0("Analysis_results/Gene_Expression/", project, "_corr_coeff_expression_all_vs_no_germline.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE)
  p <- ggscatter(df, x = "All_samples", y = "No_germline",
                 add = "reg.line",  # Add regressin line) +
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) +
    stat_cor(method = "pearson", label.x = min(df$All_samples), label.y = ymax, size = 6) +
    ggtitle(paste0(project, ": Expression", collapse = "")) +
    #xlim(-xmax,xmax) + ylim(-ymax,ymax) +
    xlab("All samples") +
    ylab("No germline") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15,face="bold"),
          axis.title.y = element_text(size = 15,face="bold"))
  
  print(p)
  dev.off()

  
}

lapply(c("BRCA", "OV", "UCEC"), test_corr)


### Correlation between fold changes with age of methylation from all samples and from samples without germline variants

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)

test_corr <- function(project){
  
  # read result from all samples
  All <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  All_coeff <- All[,c("gene", "estimate")]
  colnames(All_coeff) <- c("gene", "All_samples")
  
  # read result from samples without germline mutations
  no_germline <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age_exclude_germline.csv", collapse = ""))
  no_germline_coeff <- no_germline[,c("gene", "estimate")]
  colnames(no_germline_coeff) <- c("gene", "No_germline")
  
  # merge
  df <- merge(All_coeff, no_germline_coeff, by.x = "gene", by.y = "gene")
  
  # write source data
  write.csv(df, paste0("Source_Data/Supplementary_Fig_12a_methylation_", project, ".csv", collapse = ""), row.names = FALSE)
  
  # plot
  xmax <- max(abs(min(df$All_samples)), abs(max(df$All_samples)))
  ymax <- max(abs(min(df$No_germline)), abs(max(df$No_germline)))
  
  pdf(paste0("Analysis_results/Methylation/", project, "_corr_coeff_methylation_all_vs_no_germline.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE)
  p <- ggscatter(df, x = "All_samples", y = "No_germline",
                 add = "reg.line",  # Add regressin line) +
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) +
    stat_cor(method = "pearson", label.x = min(df$All_samples), label.y = ymax, size = 6) +
    ggtitle(paste0(project, ": Methylation", collapse = "")) +
    #xlim(-xmax,xmax) + ylim(-ymax,ymax) +
    xlab("All samples") +
    ylab("No germline") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15,face="bold"),
          axis.title.y = element_text(size = 15,face="bold"))
  
  print(p)
  dev.off()

}

lapply(c("BRCA", "OV", "UCEC"), test_corr)

