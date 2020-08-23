### Age biases in SNVs
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# keep only projects that have > 100 samples
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "LAML", "MESO", "THYM", "UCS", "UVM"))]

###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of a gene to be mutated

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# function to test association between age and mutation
test_age_mut_gene <- function(gene, project, mut_df, clinical){
  
  # get value for gene of interest
  values_gene <- mut_df[gene]
  
  values_gene$patient <- rownames(values_gene)
  rownames(values_gene) <- NULL
  colnames(values_gene) <- c("mut", "patient")
  
  # merge with age data
  df <- merge(values_gene, clinical, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- glm(mut ~ age , data = df, family = "binomial")
  summary(logit_fit)
  p.value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  coeff <- summary(logit_fit)$coefficients[,1][2]
  std.error <- summary(logit_fit)$coefficients[,2][2]
  Z <- summary(logit_fit)$coefficients[,3][2]
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  result <- result[,c("term", "estimate", "std.error", "statistic", "2.5 %", "97.5 %", "p.value")]
  colnames(result) <- c("term", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")
  result$odds <- exp(result$estimate)
  result$odds_conf.low <- exp(result$conf.low)
  result$odds_conf.high <- exp(result$conf.high)
  result <- result[,c("term", "estimate", "std.error", "conf.low", "conf.high",
                      "statistic", "odds", "odds_conf.low", "odds_conf.high", "p.value")]
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- gene
  colnames(result_df) <- c("gene", "estimate", "std.error", "conf.low", "conf.high", 
                           "statistic", "odds", "odds_conf.low", "odds_conf.high", "p.value")
  
  return(result_df)
}

# function to test mutation and age for each cancer type
age_mut <- function(project, clinical){
  
  print(paste("Start working: ", project))
  
  ### read gene-sample mutation table
  mut_df <- read.csv(paste0("Analysis_results/Mutations/1_table_mutations_samples/", project, "_mutations_filter_hypermutated.csv", collapse = ""))
  row.names(mut_df) <- mut_df$X
  mut_df$X <- NULL
  genes <- colnames(mut_df)
  
  print(paste0(project, ": ", nrow(mut_df)))  # print num samples
  
  tmp <- do.call(rbind,lapply(genes, test_age_mut_gene, project = project, mut_df = mut_df, clinical = clinical))
  rownames(tmp) <- NULL
  tmp$p.value <- as.numeric(as.character(tmp$p.value))
  tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
  tmp$Sig <- ifelse(tmp$q.value < 0.05, TRUE, FALSE)
  write.csv(tmp, paste0("Analysis_results/Mutations/2_Univariate_age_SNVs/",project, 
                        "_univariate_age_mutations_new.csv", collapse = ""), row.names = FALSE)
  return(tmp)
}

results <- lapply(projects, age_mut, clinical = clinical)
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

results_df <- results_df[,c("cancer_type", "gene", "estimate", "std.error", "conf.low", "conf.high", 
                            "statistic", "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Sig")]
write.csv(results_df, "Analysis_results/Mutations/Summary_age_SNVs_univariate_new.csv", row.names = FALSE)

