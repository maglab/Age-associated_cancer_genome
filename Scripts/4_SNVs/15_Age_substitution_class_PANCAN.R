### 6 substitution classes with age PANCAN
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(maftools)
library(ggplot2)
library(broom)
library(reshape2)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# read and extract information from maf file of each project
get_six_classes_fraction <- function(project){
  
  print(paste0("Working on: ", project))
  
  # read maf
  maf <- read.maf(maf = paste0("Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""))
  
  # remove hypermutated tumour
  #to_be_removed <- as.character(maf@variants.per.sample[maf@variants.per.sample$Variants > 1000]$Tumor_Sample_Barcode)
  
  maf@variant.classification.summary
  titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
  mut_fraction_df <- as.data.frame(titv$fraction.contribution)
  #mut_fraction_df <- mut_fraction_df[!(mut_fraction_df$Tumor_Sample_Barcode %in% to_be_removed),]
  mut_fraction_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_fraction_df$Tumor_Sample_Barcode), clean_id))
  
  # remove MSI-H tumours from COAD, READ, STAD, UCEC
  #if(project %in% c("COAD", "READ", "STAD", "UCEC")){
  #  # read MSI status file
  #  MSI <- read.csv(paste0("Data/MSI_Status/", project, "_MSI_Status.csv", collapse = ""))
  #  MSI_patients <- as.character(MSI[MSI$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)
  #  mut_fraction_df <- mut_fraction_df[!(mut_fraction_df$Tumor_Sample_Barcode %in% MSI_patients),]
  #}
  
  return(mut_fraction_df)
}

six_classes_fraction_df <- lapply(projects, get_six_classes_fraction)
six_classes_mut_fraction_all <- do.call(rbind, six_classes_fraction_df)

write.csv(six_classes_mut_fraction_all, "Analysis_results/Mutations/Mutation_6_classes.csv", row.names = FALSE)

six_classes_mut_fraction_all <- melt(six_classes_mut_fraction_all)
colnames(six_classes_mut_fraction_all) <- c("patient", "class", "fraction")

write.csv(six_classes_mut_fraction_all, "Analysis_results/Mutations/Mutation_6_classes_melted.csv", row.names = FALSE)

### linear regression for each mutation class
mut_age_class_fraction <- function(class, six_classes_mut_fraction_all, clinical){
  
  print(paste0("Working on: ", class))
  
  # get only mutation for that class
  df <- six_classes_mut_fraction_all[six_classes_mut_fraction_all$class == class,]
  
  df <- merge(df, clinical, by.x = "patient", by.y = "patient")
  
  # linear regression age and mutational burden
  lm_fit <- lm(fraction ~ age + cancer_type + gender + race, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$adj.r.squared,2)
  r_squared_1 <- summary(lm_fit)$r.squared
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  result_df <- tidy(lm_fit)
  
  write.csv(result_df, paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/PANCAN_multivariate_fraction_",class,".csv", collapse = ""), row.names = FALSE)
  
  out <- result_df[result_df$term == "age",]
  out$term <- class
  colnames(out) <- c("class", "estimate", "std.error", "statistic", "p.value")
  return(out)
}

classes <- c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")
result <- lapply(classes, mut_age_class_fraction, six_classes_mut_fraction_all, clinical)
result <- do.call(rbind, result)
result$q.value <- p.adjust(result$p.value, method = "BH")
result$Sig <- ifelse(result$q.value < 0.05, TRUE, FALSE)
write.csv(result, "Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_PANCAN_multivariate_fraction_all_classes.csv", row.names = FALSE)
