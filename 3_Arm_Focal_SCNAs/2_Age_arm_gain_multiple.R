### Age-associated arm-level recurrent gains (multiple logistic regression)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(broomExtra)
library(logistf)

### read significant arms for univariate analysis
gain_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_gain.csv")
gain_df_sig <- gain_df[gain_df$Threshold == TRUE,]

# sig cancer types
gain_projects <- as.character(unique(gain_df$cancer_type))

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# useful functions
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- gain ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- gain ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- gain ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- gain ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- gain ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- gain ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- gain ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- gain ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- gain ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- gain ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- gain ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- gain ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- gain ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- gain ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- gain ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- gain ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- gain ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- gain ~ age + purity + gender + pathologic_stage
  }
  return(model)
}


###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of an arm to be gained

# function to test association between age and arm gain
test_age_gain <- function(arm, project, purity){
  
  print(paste("Working on: ", project, ";", arm))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  ### read GISTIC arm table
  df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/broad_values_by_arm.txt", collapse = ""))
  df <- df[df$Chromosome.Arm == arm,]
  rownames(df) <- df$Chromosome.Arm
  df$Chromosome.Arm <- NULL
  colnames(df) <- unlist(lapply(colnames(df), clean_id))
  
  df <- ifelse(df > 0.25, 1, 0)   # check if gain or not
  df <- as.data.frame(t(df))
  
  df$patient <- rownames(df)
  rownames(df) <- NULL
  colnames(df) <- c("gain", "patient")
  
  # merge with age data
  df <- merge(df, clin, by.x = "patient", by.y = "patient")
  
  # merge with purity
  df <- merge(df, purity, by.x = "patient", by.y = "patient")
  
  ### select model
  model <- model_selection(project)
  
  # logistic regression
  logit_fit <- logistf(formula = model, data = df, family = "binomial")
  summary(logit_fit)
  
  result <- broomExtra::tidy_parameters(logit_fit)

  write.csv(result, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_gain/", project, "_multivariate_age_SCNA_", arm, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- arm
  result_df$cancer_type <- project
  colnames(result_df) <- c("arm", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value", "cancer_type")
  result_df <- result_df[,c("cancer_type", "arm", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")]
  
  return(result_df)
}


# function to test mutation and age for each cancer type
multivariate_result <- list()

arms <- as.character(gain_df_sig$arm)
cancer_types <- as.character(gain_df_sig$cancer_type)

for(i in 1:nrow(gain_df_sig)){
  arm <- arms[i]
  project <- cancer_types[i]
  tmp <- test_age_gain(arm = arm, project = project, purity = purity)
  multivariate_result[[i]] <- tmp
}

result_df <- do.call(rbind, multivariate_result)
result_df$q.value <- p.adjust(result_df$p.value, method = "BH")
result_df$Sig <- ifelse(result_df$q.value < 0.05, TRUE, FALSE)

write.csv(result_df, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_gain.csv", row.names = FALSE)


