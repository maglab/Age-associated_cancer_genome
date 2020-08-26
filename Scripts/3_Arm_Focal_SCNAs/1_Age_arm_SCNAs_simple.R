### Age-associated arm-level recurrent SCNAs (simple logistic regression)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(broom)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

# remove cancer type which n < 100
TCGA_projects <- TCGA_projects[!(TCGA_projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", "THYM", "UCS", "UVM"))]

###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of arm to be gained or lost

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# function to test association between age and gain of an arm
test_gain <- function(arm, values_by_arm, project){
  
  # get value for arm of interest
  values_arm <- as.data.frame(t(values_by_arm[rownames(values_by_arm) == arm,]))
  values_arm$patient <- rownames(values_arm)
  rownames(values_arm) <- NULL
  colnames(values_arm) <- c("value", "patient")
  
  # merge with age data
  df <- merge(values_arm, clinical, by.x = "patient", by.y = "patient")
  
  df$gain <- ifelse(df$value == 1, TRUE, FALSE)   # check if gain or not
  
  # logistic regression
  logit_fit <- glm(gain ~ age , data = df, family = "binomial")
  summary(logit_fit)
  p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  coeff <- summary(logit_fit)$coefficients[,1][2]
  Z <- summary(logit_fit)$coefficients[,3][2]
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  write.csv(result, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_gain/", project, "_univariate_age_recurrent_gain_", arm, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- arm
  result_df$cancer_type <- project
  colnames(result_df) <- c("arm", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "cancer_type")
  result_df <- result_df[, c("cancer_type", "arm", "estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")]
  
  return(result_df)
}

# function to test association between age and del of an arm
test_del <- function(arm, values_by_arm, project){
  
  # get value for arm of interest
  values_arm <- as.data.frame(t(values_by_arm[rownames(values_by_arm) == arm,]))
  values_arm$patient <- rownames(values_arm)
  rownames(values_arm) <- NULL
  colnames(values_arm) <- c("value", "patient")
  
  # merge with age data
  df <- merge(values_arm, clinical, by.x = "patient", by.y = "patient")
  
  df$del <- ifelse(df$value == -1, TRUE, FALSE)   # check if del or not
  
  # logistic regression
  logit_fit <- glm(del ~ age , data = df, family = "binomial")
  summary(logit_fit)
  p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  coeff <- summary(logit_fit)$coefficients[,1][2]
  Z <- summary(logit_fit)$coefficients[,3][2]
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  write.csv(result, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_del/", project, "_univariate_age_recurrent_del_", arm, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- arm
  result_df$cancer_type <- project
  colnames(result_df) <- c("arm", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "cancer_type")
  result_df <- result_df[, c("cancer_type", "arm", "estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")]
  
  return(result_df)
}

# function to test gain/del for each cancer type
age_arm <- function(project, GainOrDel){
  
  ### read significant results
  Sig_df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/broad_significance_results.txt", collapse = ""))
  colnames(Sig_df) <- c("Arm", "Num_genes", "Amp_freq", "Amp_freq_score", "Amp_z_score",
                        "Amp_q_values", "Del_freq", "Del_freq_score", "Del_z_score", "Del_q_values")
  
  sig_amp <- Sig_df[Sig_df$Amp_q_values < 0.25,]
  sig_del <- Sig_df[Sig_df$Del_q_values < 0.25,]
  
  arm_gain <- as.character(sig_amp$Arm)
  arm_del <- as.character(sig_del$Arm)
  
  ### read broad values by arm
  values_by_arm <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/broad_values_by_arm.txt", collapse = ""))
  rownames(values_by_arm) <- values_by_arm$Chromosome.Arm
  values_by_arm$Chromosome.Arm <- NULL
  colnames(values_by_arm) <- unlist(lapply(colnames(values_by_arm), clean_id))
  
  # convert values in values_by_arm by threshold (del: score < -0.25, gain: score > 0.25, no change: others)
  values_by_arm[values_by_arm < -0.25] <- -1
  values_by_arm[values_by_arm > 0.25] <- 1
  values_by_arm[values_by_arm != -1 & values_by_arm != 1] <- 0
  
  #### test association
  if(GainOrDel == "gain"){
    
    tmp <- do.call(rbind,lapply(arm_gain, test_gain, values_by_arm = values_by_arm, project = project))
    rownames(tmp) <- NULL
    tmp$p.value <- as.numeric(as.character(tmp$p.value))
    tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
    
  } else if(GainOrDel == "del"){
    
    tmp <- do.call(rbind,lapply(arm_del, test_del, values_by_arm = values_by_arm, project = project))
    rownames(tmp) <- NULL
    tmp$p.value <- as.numeric(as.character(tmp$p.value))
    tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
    
  }
  
  return(tmp)
}


### gain ##########################################################################################
all_cancers_gain_results <- lapply(TCGA_projects, age_arm, GainOrDel = "gain")
all_cancers_gain_results <- do.call(rbind, all_cancers_gain_results)

all_cancers_gain_results$Threshold <- ifelse(all_cancers_gain_results$q.value < 0.05, TRUE, FALSE)   # check if gain or not
all_cancers_gain_results$odds <- exp(all_cancers_gain_results$estimate)
all_cancers_gain_results$odds_conf.low <- exp(all_cancers_gain_results$conf.low)
all_cancers_gain_results$odds_conf.high <- exp(all_cancers_gain_results$conf.high)
all_cancers_gain_results <- all_cancers_gain_results[,c("cancer_type", "arm", "estimate", "std.error", "conf.low", "conf.high",
                                                        "statistic", "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Threshold")]
write.csv(all_cancers_gain_results, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_gain.csv", row.names = FALSE)

### del ##########################################################################################
all_cancers_del_results <- lapply(TCGA_projects, age_arm, GainOrDel = "del")
all_cancers_del_results <- do.call(rbind, all_cancers_del_results)

all_cancers_del_results$Threshold <- ifelse(all_cancers_del_results$q.value < 0.05, TRUE, FALSE)   # check if gain or not
all_cancers_del_results$odds <- exp(all_cancers_del_results$estimate)
all_cancers_del_results$odds_conf.low <- exp(all_cancers_del_results$conf.low)
all_cancers_del_results$odds_conf.high <- exp(all_cancers_del_results$conf.high)
all_cancers_del_results <- all_cancers_del_results[,c("cancer_type", "arm", "estimate", "std.error", "conf.low", "conf.high",
                                                        "statistic", "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Threshold")]

write.csv(all_cancers_del_results, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_del.csv", row.names = FALSE)

###################################################################################################
### select significant arm for further multiple logistic regression analysis
sig_gain <- all_cancers_gain_results[all_cancers_gain_results$Threshold == TRUE,]
sig_del <- all_cancers_del_results[all_cancers_del_results$Threshold == TRUE,]

write.csv(sig_gain, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_gain_sig.csv", row.names = FALSE)
write.csv(sig_del, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Univariate_age_recurrent_arm_del_sig.csv", row.names = FALSE)

