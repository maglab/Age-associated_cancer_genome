### Correlation between regression coefficient of age on methylation and regression coefficient of age on gene expression

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(ggpubr)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# remove projects that have samples < 100
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", 
                                       "READ", "THYM", "UCS", "UVM"))]
### function to perform the correlation analysis on regression coefficient of methylation with age and gene expression with age
corr_methy_exp <- function(project){
  
  ### read methylation result
  methylation <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  ### read gene expression result
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  
  genes <- as.character(methylation$gene[methylation$gene %in% expression$gene])
  length(genes)
  
  methylation <- methylation[methylation$gene %in% genes, c("gene", "estimate")]
  colnames(methylation) <- c("gene", "methylation_coeff")
  
  expression <- expression[expression$gene %in% genes, c("gene", "estimate")]
  colnames(expression) <- c("gene", "expression_coeff")
  
  df <- merge(methylation, expression, by.x = "gene", by.y = "gene")
  
  xmax <- max(abs(min(df$methylation_coeff)), abs(max(df$methylation_coeff)))
  ymax <- max(abs(min(df$expression_coeff)), abs(max(df$expression_coeff)))
  
  ### Supplementary Fig. 8
  pdf(paste0("Analysis_results/Methylation/2_Corr_coeff_methy_exp/", project, "_corr_coeff_methylation_expression.pdf", collapse = ""), width = 6, height = 4.5)
  p <- ggscatter(df, x = "methylation_coeff", y = "expression_coeff",
                 add = "reg.line",  # Add regressin line) +
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) +
    stat_cor(method = "pearson", label.x = 0.3*max(df$methylation_coeff), label.y = ymax) +
    ggtitle(project) +
    xlim(-xmax,xmax) + ylim(-ymax,ymax) +
    xlab("Regression coefficient of methylation with age") +
    ylab("Regression coefficient of gene expression with age") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"))
    
  print(p)
  dev.off()

}

lapply(projects, corr_methy_exp)

