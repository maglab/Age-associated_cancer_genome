### Age biases in mutations (multiple logistic regression)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ggplot2)
library(ggrepel)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

### results from univariate analysis
simple_regression_results <- read.csv("Analysis_results/Mutations/Summary_age_SNVs_univariate_new.csv")
simple_regression_results <- simple_regression_results[simple_regression_results$Sig == TRUE,] # only significant association between age and mutations

# cancer_types
cancer_types <- unique(as.character(simple_regression_results$cancer_type))

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- mut ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- mut ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- mut ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- mut ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- mut ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- mut ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- mut ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- mut ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- mut ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- mut ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- mut ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- mut ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- mut ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- mut ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- mut ~ age + purity + gender + pathologic_stage
  }
  return(model)
}

###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of a gene to be mutated

# function to test association between age and mutation
test_age_mut_gene <- function(gene, project, purity){
  
  print(paste("Working on: ", project, ";", gene))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  ### read gene-sample mutation table
  mut_df <- read.csv(paste0("Analysis_results/Mutations/1_table_mutations_samples/", project, "_mutations_filter_hypermutated.csv", collapse = ""))
  row.names(mut_df) <- mut_df$X
  mut_df$X <- NULL
  
  ### select model
  model <- model_selection(project)
  
  # get value for gene of interest
  values_gene <- mut_df[gene]
  
  values_gene$patient <- rownames(values_gene)
  rownames(values_gene) <- NULL
  colnames(values_gene) <- c("mut", "patient")
  
  # merge with age data
  df <- merge(values_gene, clin, by.x = "patient", by.y = "patient")
  
  # merge with purity
  df <- merge(df, purity, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- logistf(formula = model, data = df, family = "binomial")

  summary(logit_fit)
  
  result <- broomExtra::tidy_parameters(logit_fit)
  
  write.csv(result, paste0("Analysis_results/Mutations/3.1_Multivariate_age_SNVs_result_each_test/", project, "_multivariate_age_mut_", gene, "_new.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- gene
  result_df$cancer_type <- project
  colnames(result_df) <- c("gene", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value", "cancer_type")
  result_df <- result_df[,c("cancer_type", "gene", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")]
  
  return(result_df)
}


# function to test mutation and age for each cancer type
multiple_regression_result <- list()

genes <- as.character(simple_regression_results$gene)
cancer_types <- as.character(simple_regression_results$cancer_type)

for(i in 1:nrow(simple_regression_results)){
  gene <- genes[i]
  project <- cancer_types[i]
  tmp <- test_age_mut_gene(gene = gene, project = project, purity = purity)
  multiple_regression_result[[i]] <- tmp
}

result_df <- do.call(rbind, multiple_regression_result)
result_df$q.value <- p.adjust(result_df$p.value, method = "BH")
result_df$Sig <- ifelse(result_df$q.value < 0.05, TRUE, FALSE)

write.csv(result_df, "Analysis_results/Mutations/Summary_age_SNVs_multivariate_new.csv", row.names = FALSE)


### Volcano plot age-associated SNVs
df <- result_df[result_df$Sig == TRUE,]
dim(df)

# Fig. 4e
pdf("Analysis_results/Mutations/Summary_age_SNVs_multivariate_new.pdf", width = 8, height = 6) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), color = cancer_type, label = gene), data = df) +
  geom_point(size = 2.5) + 
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Cancer type-specific association between age and mutations") +
  xlim(c(-0.11,0.11)) +
  geom_label_repel(size = 3.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()


