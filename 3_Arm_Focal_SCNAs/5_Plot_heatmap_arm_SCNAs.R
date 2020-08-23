### Plot SCNA profile from GISTIC2 output
### Fig. 2f-g and Supplementary Fig. 3
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)

# read multivariate arm results
gain_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_gain.csv")
gain_df <- gain_df[gain_df$Sig == TRUE,]

del_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_del.csv")
del_df <- del_df[del_df$Sig == TRUE,]

projects <- sort(unique(c(as.character(gain_df$cancer_type), as.character(del_df$cancer_type))))

clinical <- read.csv("Data/all_clin_XML.csv")

# clean id
clean_id <- function(id){
  tmp <- unlist(strsplit(id, split = "[.]"))[1:3]
  return(paste0(tmp, collapse = "-"))
}


### Function to plot CNA profile
plot_CNA_heatmap <- function(project){
  values_by_arm <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/broad_values_by_arm.txt", collapse = ""))
  rownames(values_by_arm) <- values_by_arm$Chromosome.Arm
  values_by_arm$Chromosome.Arm <- NULL
  values_by_arm <- as.data.frame(t(values_by_arm))
  
  rownames(values_by_arm) <- unlist(lapply(rownames(values_by_arm), clean_id))
  
  values_by_arm <- values_by_arm[complete.cases(values_by_arm),]
  
  # annotation df
  annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
  annot_df <- annot_df[annot_df$patient %in% rownames(values_by_arm),]
  annot_df <- annot_df[order(annot_df$age, decreasing = TRUE),]  # sort by age
  age <- annot_df$age
  
  # annotation colours
  ha <- rowAnnotation(age = age, 
                      col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                      annotation_legend_param = list(title = "age", 
                                                     at = c(min(age), median(age), max(age)),
                                                     labels = c(min(age), median(age), max(age))))
  ha
  
  # sort df by age
  values_by_arm <- values_by_arm[as.character(annot_df$patient),]
  values_by_arm <- as.matrix(values_by_arm)
  
  values_by_arm <- values_by_arm[complete.cases(values_by_arm),]
  
  ### heatmap
  my_col = colorRamp2(c(min(values_by_arm), 0, max(values_by_arm)), c("#1d91c0", "white", "#bd0026"))
  
  pdf(paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/GISTIC_arm_heatmap/", project, "_arm_heatmap.pdf", collapse = ""), width = 8, height = 4) 
  p <- Heatmap(values_by_arm, name = "Copy Number Changes", cluster_columns = FALSE, cluster_rows = FALSE,
               show_row_names = FALSE,
               col = my_col,
               left_annotation = ha,
               column_title = project)
  print(p)
  dev.off()

}

lapply(projects, plot_CNA_heatmap)
