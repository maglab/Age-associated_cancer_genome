### Cancer-specific age-associated pathway alterations (Simple logistic regression)

library(broom)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# remove projects that have samples < 100
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", "THYM", "UCS", "UVM"))]

### pathway alteration file
pathway_df <- read.csv("Data/Pathway_alterations_clean.csv")
rownames(pathway_df) <- pathway_df$SAMPLE_BARCODE
pathway_df$SAMPLE_BARCODE <- NULL

###################################################################################################
### Logistic regression to test whether age associates with an alteration in each pathway

# function to test association between age and pathway alteration
test_age_pathway <- function(pathway, project, df_tmp, clinical){
  
  print(paste("Working on: ", project, ";", pathway))
  
  ### select data for pathway of interest
  tmp_pathway <- df_tmp[pathway]
  tmp_pathway$patient <- rownames(tmp_pathway)
  tmp_pathway <- tmp_pathway[,colnames(tmp_pathway) %in% c("patient", pathway)]
  rownames(tmp_pathway) <- NULL
  colnames(tmp_pathway) <- c("alteration", "patient")
  tmp_pathway <- tmp_pathway[complete.cases(tmp_pathway),]
  
  # merge with age data
  df <- merge(tmp_pathway, clinical, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- glm(alteration ~ age, data = df, family = "binomial")
  summary(logit_fit)
  p.value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  coeff <- summary(logit_fit)$coefficients[,1][2]
  std.error <- summary(logit_fit)$coefficients[,2][2]
  Z <- summary(logit_fit)$coefficients[,3][2]
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- pathway
  colnames(result_df) <- c("pathway", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  result_df <- result_df[, c("pathway", "estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")]
  return(result_df)
}


cancer_pathway_alterations <- function(project, pathway_df, clinical){
  
  print(paste("Working on: ", project))
  
  df_tmp <- pathway_df[rownames(pathway_df) %in% clinical[clinical$cancer_type == project,]$patient,]
  pathways <- colnames(df_tmp)
  
  tmp <- do.call(rbind,lapply(pathways, test_age_pathway, project = project, df_tmp = df_tmp, clinical = clinical))
  rownames(tmp) <- NULL
  tmp$p.value <- as.numeric(as.character(tmp$p.value))
  tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
  tmp$Sig <- ifelse(tmp$q.value < 0.05, TRUE, FALSE)
  tmp$odds <- exp(tmp$estimate)
  tmp$odds_conf.low <- exp(tmp$conf.low)
  tmp$odds_conf.high <- exp(tmp$conf.high)
  
  tmp <- tmp[,c("pathway", "estimate", "std.error", "conf.low", "conf.high", "statistic", "odds",
                "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Sig")]
  write.csv(tmp, paste0("Analysis_results/Pathway_alterations/2_Univariate_pathway_alterations/",project, 
                        "_univariate_pathway_alterations.csv", collapse = ""), row.names = FALSE)
  print(tmp)
  return(tmp)
}

results <- lapply(projects, cancer_pathway_alterations, pathway_df = pathway_df, clinical = clinical)
names(results) <- projects

results_df <- do.call(rbind, results)
results_df$cancer_type <- row.names(results_df)
row.names(results_df) <- NULL

# clean cancer type 
cancer_types <- results_df$cancer_type
clean_type <- function(cancer_type){
  return(unlist(strsplit(cancer_type, split = "[.]"))[1])
}
results_df$cancer_type <- unlist(lapply(cancer_types, clean_type))

results_df <- results_df[,c("cancer_type", "pathway", "estimate", "std.error", "conf.low", "conf.high", "statistic",
                            "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Sig")]
write.csv(results_df, "Analysis_results/Pathway_alterations/Summary_age_univariate_pathway_alterations.csv", row.names = FALSE)

