# association between SCNA score and age in each cancer type
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)
library(ggrepel)

# read data
clinical <- read.csv("Data/all_clin_XML.csv")

purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

# create a function to test the association between SCNA score and age

SCNA_age <- function(project){
  # read SCNA score
  SCNA <- read.csv(paste0("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/", project, "_normalised_SCNA_score.csv", collapse = ""))
  colnames(SCNA) <- c("patient", "rank_norm_focal", "rank_norm_chrom", "rank_norm_arm", "sum_score")
  dim(SCNA)
  df <- merge(SCNA, clinical, by.x = "patient", by.y = "patient")
  df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
  
  df$chrom_arm <- df$rank_norm_chrom + df$rank_norm_arm
  
  median_overall <- median(df$sum_score)
  median_focal <- median(df$rank_norm_focal)
  median_chromarm <- median(df$chrom_arm)
  
  print(paste0(project, " - num sample: ", nrow(df), collapse = ""))
  
  # Overall score
  lm_fit <- lm(sum_score ~ age, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_overall <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared, 2)
  coeff_overall <- summary(lm_fit)$coefficients[,1][2]
  
  result_overall <- tidy(lm_fit)
  write.csv(result_overall, paste0("Analysis_results/CNAs/2_SCNA_scores/Age_overall_SCNA_score/", project, "_univariate_age_overall_SCNA.csv", collapse = ""), row.names = FALSE)
  
  # focal score
  lm_fit <- lm(rank_norm_focal ~ age, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_focal <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared, 2)
  coeff_focal <- summary(lm_fit)$coefficients[,1][2]
  
  result_focal <- tidy(lm_fit)
  write.csv(result_focal, paste0("Analysis_results/CNAs/2_SCNA_scores/Age_focal_SCNA_score/", project, "_univariate_age_focal_SCNA.csv", collapse = ""), row.names = FALSE)
  
  # chrom arm score
  lm_fit <- lm(chrom_arm ~ age, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_chromarm <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared, 2)
  coeff_chromarm <- summary(lm_fit)$coefficients[,1][2]
  
  result_chromarm <- tidy(lm_fit)
  write.csv(result_chromarm, paste0("Analysis_results/CNAs/2_SCNA_scores/Age_chrom_arm_SCNA_score/", project, "_univariate_age_chrom_arm_SCNA.csv", collapse = ""), row.names = FALSE)
  
  result_df <- cbind(as.data.frame(result_overall[2,]), as.data.frame(result_focal[2,2:ncol(result_focal)]), as.data.frame(result_chromarm[2,2:ncol(result_chromarm)]))
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate_overall", "std.error_overall", "statistic_overall", "p.value_overall", 
                           "estimate_focal", "std.error_focal", "statistic_focal", "p.value_focal",
                           "estimate_chromarm", "std.error_chromarm", "statistic_chromarm", "p.value_chromarm")
  result_df$median_overall <- median_overall
  result_df$median_focal <- median_focal
  result_df$median_chromarm <- median_chromarm
  
  return(result_df)
}


SCNA_age_result <- lapply(TCGA_projects, SCNA_age)

SCNA_age_result_1 <- do.call(rbind, SCNA_age_result)

SCNA_age_result_1$estimate_overall <- as.numeric(as.character(SCNA_age_result_1$estimate_overall))
SCNA_age_result_1$estimate_focal <- as.numeric(as.character(SCNA_age_result_1$estimate_focal))
SCNA_age_result_1$estimate_chromarm <- as.numeric(as.character(SCNA_age_result_1$estimate_chromarm))

SCNA_age_result_1$median_overall <- as.numeric(as.character(SCNA_age_result_1$median_overall))
SCNA_age_result_1$median_focal <- as.numeric(as.character(SCNA_age_result_1$median_focal))
SCNA_age_result_1$median_chromarm <- as.numeric(as.character(SCNA_age_result_1$median_chromarm))

SCNA_age_result_1$p.value_overall <- as.numeric(as.character(SCNA_age_result_1$p.value_overall))
SCNA_age_result_1$p.value_focal <- as.numeric(as.character(SCNA_age_result_1$p.value_focal))
SCNA_age_result_1$p.value_chromarm <- as.numeric(as.character(SCNA_age_result_1$p.value_chromarm))

SCNA_age_result_1$q.value_overall <- p.adjust(SCNA_age_result_1$p.value_overall, method = "BH")
SCNA_age_result_1$q.value_focal <- p.adjust(SCNA_age_result_1$p.value_focal, method = "BH")
SCNA_age_result_1$q.value_chromarm <- p.adjust(SCNA_age_result_1$p.value_chromarm, method = "BH")

rownames(SCNA_age_result_1) <- NULL

SCNA_age_result_1$overall_threshold = as.factor(SCNA_age_result_1$q.value_overall < 0.05)
SCNA_age_result_1$focal_threshold = as.factor(SCNA_age_result_1$q.value_focal < 0.05)
SCNA_age_result_1$chromarm_threshold = as.factor(SCNA_age_result_1$q.value_chromarm < 0.05)

SCNA_age_result_2 <- SCNA_age_result_1
SCNA_age_result_2$median_overall <- NULL
SCNA_age_result_2$median_chromarm <- NULL
SCNA_age_result_2$median_focal <- NULL

write.csv(SCNA_age_result_2, "Analysis_results/CNAs/2_SCNA_scores/Summary_Univariate_Association_Age_SCNA_score.csv", row.names = FALSE)

### Multiple linear regression ##########################################################################################################################################################################
### 1) Overall scores
### Overall SCNA score ############################################################################
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- sum_score ~ age + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- sum_score ~ age + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- sum_score ~ age + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- sum_score ~ age + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- sum_score ~ age + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- sum_score ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- sum_score ~ age + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- sum_score ~ age + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- sum_score ~ age + gender + race
  } else if(project %in% c("HNSC")){
    model <- sum_score ~ age + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- sum_score ~ age + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- sum_score ~ age + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- sum_score ~ age + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- sum_score ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- sum_score ~ age + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- sum_score ~ age + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- sum_score ~ age + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- sum_score ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- sum_score ~ age + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- sum_score ~ age + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- sum_score ~ age + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- sum_score ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- sum_score ~ age + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- sum_score ~ age + gender + pathologic_stage
  }
  return(model)
}

### main function
overall_SCNAscore_age <- function(project){
  df_tmp <- read.csv(paste0("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/", project, "_normalised_SCNA_score.csv", collapse = ""))
  colnames(df_tmp) <- c("patient", "rank_norm_focal", "rank_norm_chrom", "rank_norm_arm", "sum_score")
  df_tmp$chrom_arm <- df_tmp$rank_norm_chrom + df_tmp$rank_norm_arm  # join chrom and arm scores
  df_tmp <- merge(df_tmp, purity, by.x = "patient", by.y = "patient")
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  
  model <- model_selection(project)
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$adj.r.squared, 2)
  

  ### Supplementary Fig. 2
  my_label <- paste0("adj.R-squared = ", r_squared, "\np = ", p_value)
  pdf(paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_overall_SCNA_score/", project, "_multivariate_age_overall_SCNA.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggplot(data = df_tmp, aes(x = age, y = sum_score)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    #ggtitle(paste0("Association between age and GI score in ", project, collapse = "")) +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("Overall SCNA score") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size=12,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()

  result <- tidy(lm_fit)
  
  write.csv(result, paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_overall_SCNA_score/", project, "_multivariate_age_overall_SCNA.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  
  return(result_df)
}

cancer_types <- SCNA_age_result_1[SCNA_age_result_1$q.value_overall < 0.05,]$cancer_type
results <- lapply(cancer_types, overall_SCNAscore_age)
results <- do.call(rbind, results)

# adjust p-value
results$q.value <- p.adjust(results$p.value, method = "BH")
results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)

write.csv(results, "Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_overall_SCNA_score/summary_multivariate_age_overall_SCNA.csv", row.names = FALSE)

### 2) Focal scores
### Focal SCNA score ############################################################################
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- rank_norm_focal ~ age + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- rank_norm_focal ~ age + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- rank_norm_focal ~ age + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- rank_norm_focal ~ age + gender + race
  } else if(project %in% c("HNSC")){
    model <- rank_norm_focal ~ age + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- rank_norm_focal ~ age + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- rank_norm_focal ~ age + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- rank_norm_focal ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- rank_norm_focal ~ age + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- rank_norm_focal ~ age + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- rank_norm_focal ~ age + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- rank_norm_focal ~ age + gender + pathologic_stage
  }
  return(model)
}

### main function
focal_SCNAscore_age <- function(project){
  df_tmp <- read.csv(paste0("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/", project, "_normalised_SCNA_score.csv", collapse = ""))
  colnames(df_tmp) <- c("patient", "rank_norm_focal", "rank_norm_chrom", "rank_norm_arm", "sum_score")
  df_tmp$chrom_arm <- df_tmp$rank_norm_chrom + df_tmp$rank_norm_arm  # join chrom and arm scores
  df_tmp <- merge(df_tmp, purity, by.x = "patient", by.y = "patient")
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  
  model <- model_selection(project)
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$adj.r.squared, 2)
  
  ### Supplementary Fig. 2
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  pdf(paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_focal_SCNA_score/", project, "_multivariate_age_focal_SCNA.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggplot(data = df_tmp, aes(x = age, y = rank_norm_focal)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    #ggtitle(paste0("Association between age and GI score in ", project, collapse = "")) +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("Focal SCNA score") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size=12,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result <- tidy(lm_fit)
  
  write.csv(result, paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_focal_SCNA_score/", project, "_multivariate_age_focal_SCNA.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  
  return(result_df)
}

cancer_types <- SCNA_age_result_1[SCNA_age_result_1$q.value_focal < 0.05, ]$cancer_type
results <- lapply(cancer_types, focal_SCNAscore_age)
results <- do.call(rbind, results)

# adjust p-value
results$q.value <- p.adjust(results$p.value, method = "BH")
results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)

write.csv(results, "Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_focal_SCNA_score/summary_multivariate_age_focal_SCNA.csv", row.names = FALSE)

### 3) Chrom/Arm scores
### chrom/arm SCNA score ############################################################################
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- chrom_arm ~ age + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- chrom_arm ~ age + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- chrom_arm ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- chrom_arm ~ age + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- chrom_arm ~ age + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- chrom_arm ~ age + gender + race
  } else if(project %in% c("HNSC")){
    model <- chrom_arm ~ age + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- chrom_arm ~ age + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- chrom_arm ~ age + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- chrom_arm ~ age + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- chrom_arm ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- chrom_arm ~ age + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- chrom_arm ~ age + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- chrom_arm ~ age + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- chrom_arm ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- chrom_arm ~ age + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- chrom_arm ~ age + gender + pathologic_stage
  }
  return(model)
}

### main function
chromarm_SCNAscore_age <- function(project){
  df_tmp <- read.csv(paste0("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/", project, "_normalised_SCNA_score.csv", collapse = ""))
  colnames(df_tmp) <- c("patient", "rank_norm_focal", "rank_norm_chrom", "rank_norm_arm", "sum_score")
  df_tmp$chrom_arm <- df_tmp$rank_norm_chrom + df_tmp$rank_norm_arm  # join chrom and arm scores
  df_tmp <- merge(df_tmp, purity, by.x = "patient", by.y = "patient")
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  
  model <- model_selection(project)
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$adj.r.squared, 2)
  
  ### Supplementary Fig. 2
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  pdf(paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_chrom_arm_SCNA_score/", project, "_multivariate_age_chrom_arm_SCNA.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggplot(data = df_tmp, aes(x = age, y = chrom_arm)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    #ggtitle(paste0("Association between age and GI score in ", project, collapse = "")) +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("Chrom/arm SCNA score") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size=12,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result <- tidy(lm_fit)
  
  write.csv(result, paste0("Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_chrom_arm_SCNA_score/", project, "_multivariate_age_chrom_arm_SCNA.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  
  return(result_df)
}

cancer_types <- c("KIRC", "LGG", "LUAD", "OV", "SARC", "THCA", "UCEC")
results <- lapply(cancer_types, chromarm_SCNAscore_age)
results <- do.call(rbind, results)

# adjust p-value
results$q.value <- p.adjust(results$p.value, method = "BH")
results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)

write.csv(results, "Analysis_results/CNAs/2_SCNA_scores/Multivariate_Age_chrom_arm_SCNA_score/summary_multivariate_age_chrom_arm_SCNA.csv", row.names = FALSE)


### Volcano plot ##################################################################################

### add information on whether each cancer still showed a significant association after multiple linear regression
# for overall SCNA score, sig. project included: KIRC, LGG, LUAD, OV, THCA, UCEC
SCNA_age_result_1$overall_multiple <- ifelse(SCNA_age_result_1$cancer_type %in% c("KIRC", "LGG", "LUAD", "OV", "THCA", "UCEC"), TRUE, FALSE)
# for focal SCNA score, sig. project included: KIRC, LGG, LUAD, OV, THCA, UCEC
SCNA_age_result_1$focal_multiple <- ifelse(SCNA_age_result_1$cancer_type %in% c("KIRC", "LGG", "LUAD", "OV", "THCA", "UCEC"), TRUE, FALSE)
# for chrom/arm SCNA score, sig. project included: KIRC, LGG, LUAD, OV, SARC, THCA, UCEC
SCNA_age_result_1$chromarm_multiple <- ifelse(SCNA_age_result_1$cancer_type %in% c("KIRC", "LGG", "LUAD", "OV", "SARC", "THCA", "UCEC"), TRUE, FALSE)

### add color
# overall
overall_colour <- c()
for(i in 1:nrow(SCNA_age_result_1)){
  if(SCNA_age_result_1$overall_multiple[i] == TRUE & SCNA_age_result_1$estimate_overall[i] < 0){
    overall_colour <- c(overall_colour, "Blue")
  } else if(SCNA_age_result_1$overall_threshold[i] == TRUE & SCNA_age_result_1$overall_multiple[i] == FALSE){
    overall_colour <- c(overall_colour, "Black")
  } else if(SCNA_age_result_1$overall_multiple[i] == TRUE & SCNA_age_result_1$estimate_overall[i] > 0){
    overall_colour <- c(overall_colour, "Red")
  } else {
    overall_colour <- c(overall_colour, "Grey")
  }
}

SCNA_age_result_1$overall_colour <- factor(overall_colour)

# focal
focal_colour <- c()
for(i in 1:nrow(SCNA_age_result_1)){
  if(SCNA_age_result_1$focal_multiple[i] == TRUE & SCNA_age_result_1$estimate_focal[i] < 0){
    focal_colour <- c(focal_colour, "Blue")
  } else if(SCNA_age_result_1$focal_threshold[i] == TRUE & SCNA_age_result_1$focal_multiple[i] == FALSE){
    focal_colour <- c(focal_colour, "Black")
  } else if(SCNA_age_result_1$focal_multiple[i] == TRUE & SCNA_age_result_1$estimate_focal[i] > 0){
    focal_colour <- c(focal_colour, "Red")
  } else {
    focal_colour <- c(focal_colour, "Grey")
  }
}

SCNA_age_result_1$focal_colour <- factor(focal_colour)

# chrom/arm
chromarm_colour <- c()
for(i in 1:nrow(SCNA_age_result_1)){
  if(SCNA_age_result_1$chromarm_multiple[i] == TRUE & SCNA_age_result_1$estimate_chromarm[i] < 0){
    chromarm_colour <- c(chromarm_colour, "Blue")
  } else if(SCNA_age_result_1$chromarm_threshold[i] == TRUE & SCNA_age_result_1$chromarm_multiple[i] == FALSE){
    chromarm_colour <- c(chromarm_colour, "Black")
  } else if(SCNA_age_result_1$chromarm_multiple[i] == TRUE & SCNA_age_result_1$estimate_chromarm[i] > 0){
    chromarm_colour <- c(chromarm_colour, "Red")
  } else {
    chromarm_colour <- c(chromarm_colour, "Grey")
  }
}

SCNA_age_result_1$chromarm_colour <- factor(chromarm_colour)

### Fig. 2a (Overall)
pdf("Analysis_results/CNAs/2_SCNA_scores/Summary_multivariate_age_overall_SCNA.pdf", width = 6, height = 4.5) 
p <- ggplot(aes(x = estimate_overall, y = -log10(q.value_overall), size = median_overall, 
                color = overall_colour, label = cancer_type), data = SCNA_age_result_1) +
  geom_point() +
  scale_color_manual(values = c("#000000", "#1d91c0", "#bdbdbd", "#a50f15")) +
  scale_size_continuous(range = c(1, 4)) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Overall SCNA score") +
  xlim(c(-0.025,0.025)) +
  geom_label_repel(size = 3) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()


### Fig. 2b (Chrom/Arm level)
pdf("Analysis_results/CNAs/2_SCNA_scores/Summary_multivariate_age_chromarm_SCNA.pdf", width = 6, height = 4.5) 
p <- ggplot(aes(x = estimate_chromarm, y = -log10(q.value_chromarm), size = median_chromarm, 
                color = chromarm_colour, label = cancer_type), data = SCNA_age_result_1) +
  geom_point() +
  scale_color_manual(values = c("#1d91c0", "#bdbdbd", "#a50f15")) +
  scale_size_continuous(range = c(1, 4)) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Chromosome/arm SCNA score") +
  xlim(c(-0.015,0.015)) +
  geom_label_repel(size = 3) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()

### Fig. 2c (Focal level)
pdf("Analysis_results/CNAs/2_SCNA_scores/Summary_multivariate_age_focal_SCNA.pdf", width = 6, height = 4.5) 
p <- ggplot(aes(x = estimate_focal, y = -log10(q.value_focal), size = median_focal, 
                color = focal_colour, label = cancer_type), data = SCNA_age_result_1) +
  geom_point() +
  scale_color_manual(values = c("#000000", "#1d91c0", "#bdbdbd", "#a50f15")) +
  scale_size_continuous(range = c(1, 4)) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Focal SCNA score") +
  xlim(c(-0.01,0.01)) +
  geom_label_repel(size = 3) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()

# here in the figures are the plot between 
# - X axis: simple linear regression coefficient
# - Y axis: simple linear regression -log10(adjusted p value by BH method)
# Colour represents significant (multiple linear adjusted p value < 0.05)
# - Grey: not significant
# - Red: significant and positive coefficient
# - Blue: significant and negative coefficient
# - Black: significant in univariate but not in multivariate analyses
# Size corresponses to median SCNA score.

