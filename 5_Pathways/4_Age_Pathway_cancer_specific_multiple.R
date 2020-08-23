### Cancer-specific age-associated pathway alterations (Multiple logistic regression)

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(ggplot2)
library(reshape2)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

### results from simple logistic regression analysis
simple_regression_results <- read.csv("Analysis_results/Pathway_alterations/Summary_age_univariate_pathway_alterations.csv")
simple_regression_results <- simple_regression_results[simple_regression_results$Sig == TRUE,] # only significant association between age and mutations
dim(simple_regression_results)   # 35 alterations

# cancer_types
cancer_types <- unique(as.character(simple_regression_results$cancer_type))  # 16 cancer types

### pathway alteration file
pathway_df <- read.csv("Data/Pathway_alterations_clean.csv")
rownames(pathway_df) <- pathway_df$SAMPLE_BARCODE
pathway_df$SAMPLE_BARCODE <- NULL

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- alteration ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- alteration ~ age + purity + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- alteration ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- alteration ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- alteration ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- alteration ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- alteration ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- alteration ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- alteration ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- alteration ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- alteration ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- alteration ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- alteration ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- alteration ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- alteration ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- alteration ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- alteration ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- alteration ~ age + purity + gender + pathologic_stage
  }
  return(model)
}


###################################################################################################
### Logistic regression to test whether age associates with an alteration in each pathway

# function to test association between age and pathway alteration
test_age_pathway <- function(pathway, project, pathway_df, purity){
  
  print(paste("Working on: ", project, ";", pathway))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  ### select data for pathway of interest
  tmp_pathway <- pathway_df[as.character(rownames(pathway_df)) %in% as.character(clin$patient), ]
  tmp_pathway$patient <- rownames(tmp_pathway)
  tmp_pathway <- tmp_pathway[,colnames(tmp_pathway) %in% c("patient", pathway)]
  rownames(tmp_pathway) <- NULL
  colnames(tmp_pathway) <- c("alteration", "patient")
  tmp_pathway <- tmp_pathway[complete.cases(tmp_pathway),]
  
  ### select model
  model <- model_selection(project)
  
  # merge with age data
  df <- merge(tmp_pathway, clin, by.x = "patient", by.y = "patient")
  
  # merge with purity
  df <- merge(df, purity, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- glm(formula = model, data = df, family = "binomial")
  
  summary(logit_fit)
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  result <- result[,c("term", "estimate", "std.error", "statistic", "2.5 %", "97.5 %", "p.value")]
  colnames(result) <- c("term", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")
  write.csv(result, paste0("Analysis_results/Pathway_alterations/3.1_Multivariate_pathway_alterations_each_test/", project, "_pathway_alterations_multivariate_", pathway, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- pathway
  colnames(result_df) <- c("pathway", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")
  
  result_df$cancer_type <- project
  result_df <- result_df[,c("cancer_type", "pathway", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")]
  
  return(result_df)
}

# test pathway alterations and age for each cancer type
multiple_regression_result <- list()

pathways <- as.character(simple_regression_results$pathway)
cancer_types <- as.character(simple_regression_results$cancer_type)

for(i in 1:nrow(simple_regression_results)){
  pathway <- pathways[i]
  project <- cancer_types[i]
  tmp <- test_age_pathway(pathway = pathway, project = project, pathway_df = pathway_df, purity = purity)
  multiple_regression_result[[i]] <- tmp
}

result_df <- do.call(rbind, multiple_regression_result)
result_df$q.value <- p.adjust(result_df$p.value, method = "BH")
result_df$Sig <- ifelse(result_df$q.value < 0.05, TRUE, FALSE)

result_df$odds <- exp(result_df$estimate)
result_df$odds_conf.low <- exp(result_df$conf.low)
result_df$odds_conf.high <- exp(result_df$conf.high)

result_df <- result_df[,c("cancer_type", "pathway", "estimate", "std.error", "conf.low", "conf.high", "statistic",
                          "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Sig")]
write.csv(result_df, "Analysis_results/Pathway_alterations/Summary_age_multivariate_pathway_alterations.csv", row.names = FALSE)

### Fig. 5b
df <- result_df[result_df$Sig == TRUE,]

df$direction <- ifelse(df$estimate > 0, "increase", "decrease")
cols <- c("decrease" = "#1D91C0", "increase" = "#a50f15")

pdf("Analysis_results/Pathway_alterations/Summary_age_multivariate_pathway_alterations_dot_plot.pdf", width = 8, height = 4) 
p <- ggplot(data = df, aes(x = pathway, y = cancer_type)) +
  geom_point(aes(fill = direction, size = -log10(q.value)), colour="black", pch=21) +
  scale_x_discrete(limits = names(sort(table(as.character(df$pathway)), decreasing = TRUE))) +
  scale_y_discrete(limits = names(sort(table(as.character(df$cancer_type)), decreasing = FALSE))) +
  scale_fill_manual(values = cols) +
  scale_size_continuous(range = c(2, 5)) +
  ggtitle("Cancer-specific association between age and pathway alterations") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9")) +
  guides(size = guide_legend("-log10(adj. p-val)", order = 1))
print(p)
dev.off()


