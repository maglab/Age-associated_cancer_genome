### plot heatmap of focal-level CNAs of the regions that significantly associated with age
### Fig. 3b-c and Supplementary Fig. 4
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)

# read significant results
gain <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain_sig.csv")

loss <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_loss_sig.csv")

cancer_types <- c(as.character(gain$cancer_type), as.character(loss$cancer_type))
cancer_types <- unique(cancer_types)

# clean id
clean_id <- function(id){
  tmp <- unlist(strsplit(id, split = "[.]"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# function to plot focal-level CNAs status by age
plot_CNAs <- function(project){
  
  print(paste0("Start working on: ", project))
  
  # read clinical files
  clinical <- read.csv("Data/all_clin_XML.csv")
  
  # cancer of interest data
  gain_df <- gain[gain$cancer_type == project,]
  gain_df$GainOrLoss <- rep("gain", nrow(gain_df))
  gain_df$direction <- ifelse(gain_df$estimate > 0, "increase", "decrease")
  
  loss_df <- loss[loss$cancer_type == project,]
  loss_df$GainOrLoss <- rep("loss", nrow(loss_df))
  loss_df$direction <- ifelse(loss_df$estimate > 0, "increase", "decrease")
  
  if(nrow(gain_df) > 0 & nrow(loss_df) > 0){
    my_df <- rbind(gain_df, loss_df)
  } else if(nrow(gain_df) > 0){
    my_df <- gain_df
  } else {
    my_df <- loss_df
  }
  
  print(my_df)
  
  # get all peaks that are significant associated with age
  gain_peaks <- as.character(gain[gain$cancer_type == project,]$Unique.Name)
  loss_peaks <- as.character(loss[loss$cancer_type == project,]$Unique.Name)
  sig_peaks <- c(gain_peaks, loss_peaks)
  
  # read GISTIC result file for the project of interest
  values_focal <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
  values_focal$X <- NULL   # remove last column which has nothing but NAs
  
  values_focal <- values_focal[values_focal$Unique.Name %in% sig_peaks,]  # keep only age-associated regions
  values_focal <- values_focal[, -(3:9)]
  
  # clean ID
  colnames(values_focal) <- c(colnames(values_focal)[1:2], unlist(lapply(colnames(values_focal)[3:ncol(values_focal)], clean_id)))
  
  rownames(values_focal) <- values_focal$Descriptor
  values_focal$Descriptor <- NULL
  values_focal$Unique.Name <- NULL
  
  # annotation df
  annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
  annot_df <- annot_df[annot_df$patient %in% colnames(values_focal),]
  annot_df <- annot_df[order(annot_df$age, decreasing = FALSE),]  # sort by age
  age <- annot_df$age
  
  # annotation colours
  ha <- columnAnnotation(age = age, 
                         col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                         annotation_legend_param = list(title = "age", 
                                                        at = c(min(age), median(age), max(age)),
                                                        labels = c(min(age), median(age), max(age))),
                         annotation_name_side = "left")
  ha
  
  # annotation increase or decrease freq with age
  annot_direction <- my_df[,c("region", "direction", "GainOrLoss"),]
  
  row_ha <- rowAnnotation(direction = annot_direction$direction,
                          gain_or_loss = annot_direction$GainOrLoss,
                          col = list(direction = c("increase" = "#fdae61", "decrease" = "#7bccc4"),
                                     gain_or_loss = c("gain" = "#a50026", "loss" = "#313695")))
  
  # sort df by age
  values_focal <- values_focal[,as.character(annot_df$patient)]
  values_focal <- as.matrix(values_focal)
  
  #values_focal <- values_focal[complete.cases(values_focal),]
  
  ### heatmap
  my_col = colorRamp2(c(min(values_focal), 0, max(values_focal)), c("#1d91c0", "white", "#bd0026"))
  
  pdf(paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/GISTIC_focal_heatmap/", project, "_focal_heatmap.pdf", collapse = ""), width = 10, height = 6)
  p <- Heatmap(values_focal, name = "Copy Number Changes", cluster_columns = FALSE, cluster_rows = FALSE,
               show_row_names = TRUE,
               row_names_side = "left",
               show_column_names = FALSE,
               col = my_col,
               bottom_annotation = ha,
               right_annotation = row_ha,
               column_title = project)
  
  print(p)
  dev.off()

  
}

lapply(cancer_types[cancer_types != "UCEC"], plot_CNAs)   # UCEC has the region that sig. for both gain and loss

############################################## UCEC ###############################################
project <- "UCEC"

# read clinical files
clinical <- read.csv("Data/all_clin_XML.csv")

# cancer of interest data
gain_df <- gain[gain$cancer_type == project,]
gain_df$GainOrLoss <- rep("gain", nrow(gain_df))
gain_df$direction <- ifelse(gain_df$estimate > 0, "increase", "decrease")

loss_df <- loss[loss$cancer_type == project,]
loss_df$GainOrLoss <- rep("loss", nrow(loss_df))
loss_df$direction <- ifelse(loss_df$estimate > 0, "increase", "decrease")

# get all peaks that are significant associated with age
gain_peaks <- as.character(gain[gain$cancer_type == project,]$Unique.Name)
loss_peaks <- as.character(loss[loss$cancer_type == project,]$Unique.Name)

### Gain ##########################################################################################
# read GISTIC result file for the project of interest
values_focal <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
values_focal$X <- NULL   # remove last column which has nothing but NAs

values_focal <- values_focal[values_focal$Unique.Name %in% gain_peaks,]  # keep only age-associated regions
values_focal <- values_focal[, -(3:9)]

# clean ID
colnames(values_focal) <- c(colnames(values_focal)[1:2], unlist(lapply(colnames(values_focal)[3:ncol(values_focal)], clean_id)))

rownames(values_focal) <- values_focal$Descriptor
values_focal$Descriptor <- NULL
values_focal$Unique.Name <- NULL

# annotation df
annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
annot_df <- annot_df[annot_df$patient %in% colnames(values_focal),]
annot_df <- annot_df[order(annot_df$age, decreasing = FALSE),]  # sort by age
age <- annot_df$age

# annotation colours
ha <- columnAnnotation(age = age, 
                       col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                       annotation_legend_param = list(title = "age", 
                                                      at = c(min(age), median(age), max(age)),
                                                      labels = c(min(age), median(age), max(age))),
                       annotation_name_side = "left")
ha

# annotation increase or decrease freq with age
annot_direction <- gain_df[,c("region", "direction", "GainOrLoss"),]

row_ha <- rowAnnotation(direction = annot_direction$direction,
                        gain_or_loss = annot_direction$GainOrLoss,
                        col = list(direction = c("increase" = "#fdae61", "decrease" = "#7bccc4"),
                                   gain_or_loss = c("gain" = "#a50026", "loss" = "#313695")))

# sort df by age
values_focal <- values_focal[,as.character(annot_df$patient)]
values_focal <- as.matrix(values_focal)

#values_focal <- values_focal[complete.cases(values_focal),]

### heatmap
my_col = colorRamp2(c(min(values_focal), 0, max(values_focal)), c("#1d91c0", "white", "#bd0026"))

pdf(paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/GISTIC_focal_heatmap/", project, "_focal_heatmap_gain.pdf", collapse = ""), width = 10, height = 6)
p <- Heatmap(values_focal, name = "Copy Number Changes", cluster_columns = FALSE, cluster_rows = FALSE,
             show_row_names = TRUE,
             row_names_side = "left",
             show_column_names = FALSE,
             col = my_col,
             bottom_annotation = ha,
             right_annotation = row_ha,
             column_title = "UCEC Gain")

print(p)
dev.off()


### Loss ##########################################################################################
# read GISTIC result file for the project of interest
values_focal <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
values_focal$X <- NULL   # remove last column which has nothing but NAs

values_focal <- values_focal[values_focal$Unique.Name %in% loss_peaks,]  # keep only age-associated regions
values_focal <- values_focal[, -(3:9)]

# clean ID
colnames(values_focal) <- c(colnames(values_focal)[1:2], unlist(lapply(colnames(values_focal)[3:ncol(values_focal)], clean_id)))

rownames(values_focal) <- values_focal$Descriptor
values_focal$Descriptor <- NULL
values_focal$Unique.Name <- NULL

# annotation df
annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
annot_df <- annot_df[annot_df$patient %in% colnames(values_focal),]
annot_df <- annot_df[order(annot_df$age, decreasing = FALSE),]  # sort by age
age <- annot_df$age

# annotation colours
ha <- columnAnnotation(age = age, 
                       col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                       annotation_legend_param = list(title = "age", 
                                                      at = c(min(age), median(age), max(age)),
                                                      labels = c(min(age), median(age), max(age))),
                       annotation_name_side = "left")
ha

# annotation increase or decrease freq with age
annot_direction <- loss_df[,c("region", "direction", "GainOrLoss"),]

row_ha <- rowAnnotation(direction = annot_direction$direction,
                        gain_or_loss = annot_direction$GainOrLoss,
                        col = list(direction = c("increase" = "#fdae61", "decrease" = "#7bccc4"),
                                   gain_or_loss = c("gain" = "#a50026", "loss" = "#313695")))

# sort df by age
values_focal <- values_focal[,as.character(annot_df$patient)]
values_focal <- as.matrix(values_focal)

#values_focal <- values_focal[complete.cases(values_focal),]

### heatmap
my_col = colorRamp2(c(min(values_focal), 0, max(values_focal)), c("#1d91c0", "white", "#bd0026"))

pdf(paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/GISTIC_focal_heatmap/", project, "_focal_heatmap_loss.pdf", collapse = ""), width = 10, height = 6)
p <- Heatmap(values_focal, name = "Copy Number Changes", cluster_columns = FALSE, cluster_rows = FALSE,
             show_row_names = TRUE,
             row_names_side = "left",
             show_column_names = FALSE,
             col = my_col,
             bottom_annotation = ha,
             right_annotation = row_ha,
             column_title = "UCEC Loss")

print(p)
dev.off()

