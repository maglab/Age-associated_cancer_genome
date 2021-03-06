### Mutational burden with age (cancer type-specific)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(ggplot2)
library(broom)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to read mutational burden file and plot with age
mut_burden <- function(project){
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id))
  mut_burden_df <- merge(mut_burden_df, clinical, by.x = "Tumor_Sample_Barcode", by.y = "patient")
  mut_burden_df$log_mut <- log10(mut_burden_df$total)
  
  print(nrow(mut_burden_df))
  
  # linear regression age and mutational burden
  lm_fit <- lm(log_mut ~ age, data=mut_burden_df) 
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared,2)
  r_squared_1 <- summary(lm_fit)$r.squared
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  result_df <- tidy(lm_fit)
  result_df <- as.data.frame(result_df[result_df$term == "age",])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  result_df$Rsquared <- r_squared_1
  result_df <- result_df[, c("cancer_type", "Rsquared", "estimate", "std.error", "statistic", "p.value")]
  colnames(result_df) <- c("cancer_type", "R-squared", "estimate", "std.error", "statistic", "p.value")
  return(result_df)
}

mut_burden_result <- lapply(projects, mut_burden)
mut_burden_result <- do.call(rbind, mut_burden_result)

mut_burden_result$q.value <- p.adjust(mut_burden_result$p.value, method = "BH")
mut_burden_result$Sig <- ifelse(mut_burden_result$q.value < 0.05, TRUE, FALSE)

write.csv(mut_burden_result, "Analysis_results/Mutations/Summary_mut_burden_with_age.csv", row.names = FALSE)


### Multiple linear regression ############################################################################################################################

projects <- as.character(mut_burden_result[mut_burden_result$Sig == TRUE, ]$cancer_type)

# sub-function to select a regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- log_mut ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- log_mut ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- log_mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- log_mut ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- log_mut ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- log_mut ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- log_mut ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- log_mut ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- log_mut ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- log_mut ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- log_mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- log_mut ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- log_mut ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- log_mut ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- log_mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- log_mut ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- log_mut ~ age + purity + gender + pathologic_stage
  }
  return(model)
}

### function to read mutational burden file and plot with age
mut_burden_multiple <- function(project, purity){
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df <- mut_burden_df[,c("Tumor_Sample_Barcode", "total")]
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id))
  mut_burden_df$log_mut <- log10(mut_burden_df$total)
  print(nrow(mut_burden_df))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  mut_burden_df <- merge(mut_burden_df, clin, by.x = "Tumor_Sample_Barcode", by.y = "patient")
  
  mut_burden_df <- merge(mut_burden_df, purity, by.x = "Tumor_Sample_Barcode", by.y = "patient")
  
  # write source data
  write.csv(mut_burden_df, paste0("Source_Data/Supplementary_Fig_7a_", project, ".csv", collapse = ""), row.names = FALSE)
  
  model <- model_selection(project)
  
  # linear regression age and mutational burden
  lm_fit <- lm(formula = model, data=mut_burden_df) 
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$r.squared,2)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$adj.r.squared,2)
  r_squared_1 <- summary(lm_fit)$adj.r.squared
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  ### Supplementary Fig. 5b
  pdf(paste0("Analysis_results/Mutations/4_Mut_burden_plot/", project, "_mut_burden_multivariate.pdf"), width = 6, height = 4.5, useDingbats=FALSE) 
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  p <- ggplot(data = mut_burden_df, aes(x = age, y = log_mut)) + 
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("log10(total mutations)") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face="bold"),
          panel.background = element_blank(),
          #panel.border = element_rect(linetype = "solid", fill = NA),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 6,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result_df <- tidy(lm_fit)
  write.csv(result_df, paste0("Analysis_results/Mutations/4_Mut_burden_plot/", project, "_mut_burden_multivariate.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result_df[result_df$term == "age",])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  result_df$Rsquared <- r_squared_1
  result_df <- result_df[, c("cancer_type", "Rsquared", "estimate", "std.error", "statistic", "p.value")]
  colnames(result_df) <- c("cancer_type", "adj. R-squared", "estimate", "std.error", "statistic", "p.value")
  return(result_df)
}

df <- lapply(projects, mut_burden_multiple, purity = purity)
df <- do.call(rbind, df)
df$q.value <- p.adjust(df$p.value, method = "BH")

df$Sig <- ifelse(df$q.value < 0.05, TRUE, FALSE)
write.csv(df, "Analysis_results/Mutations/Summary_mut_burden_with_age_multivariate.csv", row.names = FALSE)


########## UCEC exclude MSI-H and POLE/POLD mutations ###########################################################################
# mutational burden
mut_burden <- read.csv("Analysis_results/Mutations/Mut_burden/UCEC_mutational_burdens.csv")
mut_burden <- mut_burden[,c("Tumor_Sample_Barcode", "total")]
mut_burden$Tumor_Sample_Barcode <- unlist(lapply(mut_burden$Tumor_Sample_Barcode,clean_id))

# clinical
clin <- read.csv("Data/clinical_XML_interest/TCGA-UCEC_clinical_XML_interest.csv")

# merge clin and mut_burden
df <- merge(clin, mut_burden, by.x = "patient", by.y = "Tumor_Sample_Barcode")
df <- merge(df, purity, by.x = "patient", by.y = "patient")

# MSI status
MSI <- read.csv("Data/MSI_Status/UCEC_MSI_Status.csv")
MSI <- MSI[,c("bcr_patient_barcode", "mononucleotide_and_dinucleotide_marker_panel_analysis_status")]
colnames(MSI) <- c("patient", "msi")

# consider only patients with MSI information
shared_samples <- as.character(df$patient[df$patient %in% MSI$patient])
MSI <- MSI[MSI$patient %in% shared_samples,]
MSI_H_patients <- as.character(MSI[MSI$msi == "MSI-H",]$patient)

df <- df[df$patient %in% shared_samples,]

# PLOE/POLD1 mutations
POLE_POLD1 <- read.csv("Analysis_results/Mutations/6_POLE_POLD1/UCEC_pol_mutations.csv")
POLE_POLD1_patients <- as.character(POLE_POLD1[POLE_POLD1$POLD1 == 1 | POLE_POLD1$POLE == 1,]$X)


### samples to be removed
samples_to_be_removed <- unique(c(MSI_H_patients, POLE_POLD1_patients))

df_removed <- df[!(df$patient %in% samples_to_be_removed),]
dim(df_removed)

### linear regression
df_removed$log_mut <- log10(df_removed$total)

# write source data
write.csv(df_removed, "Source_Data/Supplementary_Fig_7c.csv", row.names = FALSE)

# linear regression age and mutational burden: multivariate
lm_fit <- lm(log_mut ~ age + purity + race + figo_stage + histologic_grade, data=df_removed) 
p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
r_squared <- round(summary(lm_fit)$adj.r.squared,2)
r_squared_1 <- summary(lm_fit)$adj.r.squared
coeff <- summary(lm_fit)$coefficients[,1][2]

pdf("Analysis_results/Mutations/UCEC_mut_burden_no_MSI-H_POLE_POLD1_multivariate.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
p <- ggplot(data = df_removed, aes(x = age, y = log_mut)) + 
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE) +
  ggtitle("UCEC (excluded MSI-H, POLE and POLD1 mutations)") +
  xlab("Age at diagnosis") +
  ylab("log10(total mutations)") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face="bold"),
        panel.background = element_blank(),
        #panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 6,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()


