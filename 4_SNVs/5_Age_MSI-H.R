### age-associated MSI status

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ggplot2)

projects <- c("UCEC", "COAD", "READ", "STAD")

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# sub-function to select a regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("UCEC")){
    model <- MSI ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("COAD", "READ")){
    model <- MSI ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("STAD")){
    model <- MSI ~ age + purity + gender + pathologic_stage + histologic_grade
  } 
  return(model)
}

### main function
MSI_H_age <- function(project){
  ### mut burden
  mut_burden <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden$Tumor_Sample_Barcode), clean_id))
  mut_burden <- mut_burden[,c("Tumor_Sample_Barcode", "total")]
  
  ### MSI status file
  MSI <- read.csv(paste0("Data/MSI_Status/", project, "_MSI_Status.csv", collapse = ""))
  MSI <- MSI[,c("bcr_patient_barcode", "mononucleotide_and_dinucleotide_marker_panel_analysis_status")]
  colnames(MSI) <- c("patient", "msi")
  MSI <- MSI[MSI$patient %in% mut_burden$Tumor_Sample_Barcode,]
  
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  
  df <- merge(clin, MSI, by.x = "patient", by.y = "patient")
  df <- merge(df, mut_burden, by.x = "patient", by.y = "Tumor_Sample_Barcode")
  head(df)
  
  df <- merge(df, purity, by.x = "patient", by.y = "patient")
  
  df$MSI <- ifelse(df$msi == "MSI-H", TRUE, FALSE)
  
  model <- model_selection(project)
  
  # logistic regression
  logit_fit <- logistf(formula = model, data = df, family = "binomial")
  summary(logit_fit)
  result <- broomExtra::tidy_parameters(logit_fit)
  p_value <- round(result$p.value[2],4)
  write.csv(result, "Analysis_results/Mutations/5_UCEC_low_burden_age_SNVs/UCEC_MSI_with_age.csv", row.names = FALSE)
  
  ### Fig. 4b and Supplementary Fig. 6a
  my_label <- paste0("p = ", p_value)
  pdf("Analysis_results/Mutations/5_UCEC_low_burden_age_SNVs/UCEC_MSI_with_age.pdf", width = 3, height = 4) 
  p <- ggplot(aes(x = MSI, y = age, fill = MSI), data = df) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(width = 0.4, fill = "white") +
    ggtitle("UCEC: MSI-H and age") +
    xlab("MSI-H") +
    ylab("Age at diagnosis") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size=12,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          legend.title = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          #panel.border = element_rect(linetype = "solid", fill = NA),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
}

lapply(projects, MSI_H_age)
