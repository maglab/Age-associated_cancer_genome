### Age biases in focal-level recurrent SCNAs (multiple logistic regression)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(broomExtra)
library(logistf)

### read significant focal regions for simple logistic regression analysis
focal_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Univariate_age_recurrent_focal_sig.csv")
gain_df <- focal_df[focal_df$GainOrLoss == "gain",]

# sig cancer types
gain_projects <- sort(as.character(unique(gain_df$cancer_type)))

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
### Logistic regression to test whether age associates with increased/decreased possibility of a focal region to be gained

# function to test association between age and focal region gain
test_age_gain <- function(Unique.Name, project, purity){
  
  print(paste("Working on: ", project, ";", Unique.Name))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  ### read GISTIC region table
  df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
  df$X <- NULL
  df$Descriptor <- unlist(lapply(as.character(df$Descriptor), gsub, pattern = " ", replacement = "", fixed = TRUE))
  df <- df[df$Unique.Name == Unique.Name,]
  
  # matching df for peak name and regions
  match_df <- df[,c("Unique.Name", "Descriptor")]
  rownames(match_df) <- NULL                 
  
  df <- df[, -(7:9)]    # remove some columns
  df <- df[,-(3:6)]   # remove some columns to keep only columns for samples
  
  rownames(df) <- df$Unique.Name
  df$Unique.Name <- NULL
  df$Descriptor <- NULL
  
  colnames(df) <- unlist(lapply(colnames(df), clean_id))
  
  ### convert dataframe to 1, 0 (TRUE, FALSE)
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
  
  Descriptor <- match_df$Descriptor
  Name <- unlist(strsplit(Unique.Name, split = " "))
  Name <- Name[Name != ""]
  Name <- paste0(Name[c(1,2,3)], collapse = "_")
  Name <- paste0(Name, "_", Descriptor, collapse = "")
  write.csv(result, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain/", project, "_multivariate_age_SCNA_", Name, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- Descriptor
  result_df$cancer_type <- project
  result_df$Unique.Name <- Unique.Name
  colnames(result_df) <- c("region", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value", "cancer_type", "Unique.Name")
  result_df <- result_df[,c("cancer_type", "Unique.Name", "region", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")]
  
  return(result_df)
}


# function to test focal SCNA and age for each cancer type
multivariate_result <- list()

regions <- as.character(gain_df$Unique.Name)
cancer_types <- as.character(gain_df$cancer_type)

for(i in 1:nrow(gain_df)){
  Unique.Name <- regions[i]
  project <- cancer_types[i]
  tmp <- test_age_gain(Unique.Name = Unique.Name, project = project, purity = purity)
  multivariate_result[[i]] <- tmp
}

result_df <- do.call(rbind, multivariate_result)
result_df$q.value <- p.adjust(result_df$p.value, method = "BH")
result_df$Sig <- ifelse(result_df$q.value < 0.05, TRUE, FALSE)

write.csv(result_df, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain.csv", row.names = FALSE)

# select only sig regions
result_df <- result_df[result_df$Sig == TRUE,]
# remove duplicated peaks by keep the most significant one
duplicated_peaks <- result_df[which(duplicated(result_df[,c("cancer_type", "region")]) == TRUE),]   # and UCEC 3q26.2
result_df <- result_df[!(result_df$cancer_type == "UCEC" & as.character(result_df$Unique.Name) == "Amplification Peak 11 - CN values"),]

write.csv(result_df, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain_sig.csv", row.names = FALSE)

