# association between WGD status and age in each cancer type
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)

# read data
samples_of_interest <- read.csv("Data/samples_in_ASCAT_and_seg.csv")

GI <- read.table("Data/TCGA.giScores.wgd.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(GI, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df)

cancer_types <- sort(as.character(unique(df$cancer_type)))

# create a function to test the association between WGD status and age
WGD_age <- function(cancer_type){
  df <- df[df$cancer_type == cancer_type,]
  print(paste0(cancer_type, " - num sample: ", nrow(df), collapse = ""))
  logit_fit <- glm(wgd ~ age , data = df, family = "binomial")
  p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(logit_fit)$coefficients[,4][2])
  coeff <- summary(logit_fit)$coefficients[,1][2]
  
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
  
  write.csv(result, paste0("Analysis_results/Structural_Alterations/3_Age_WGD/", cancer_type, "_univariate_age_WGD.csv", collapse = ""), row.names = FALSE)
  
  pdf(paste0("Analysis_results/Structural_Alterations/3_Age_WGD/", cancer_type, "_univariate_age_WGD.pdf", collapse = ""), width = 3, height = 4, useDingbats=FALSE) 
  my_label <- paste0("p = ", p_value)
  p <- ggplot(df, aes(x=wgd, y=age, fill=wgd)) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(width = 0.4, fill = "white") +
    ggtitle(cancer_type) +
    xlab("Whole Genome Duplication") +
    ylab("Age at diagnosis") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size=14,face="bold"),
          axis.title.y = element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()

  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- cancer_type
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "conf.low", "conf.high", "statistic",
                           "odds", "odds_conf.low", "odds_conf.high", "p.value")
  return(result_df)
}

results <- lapply(cancer_types, WGD_age)
results <- do.call(rbind, results)
results$q.value <- p.adjust(results$p.value, method = "BH")
results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)

write.csv(results, "Analysis_results/Structural_Alterations/3_Age_WGD/Summary_univariate_association_age_WGD.csv", row.names = FALSE)

for_multiple_regression <- results[results$Sig == TRUE,]$cancer_type

### Multiple logistic regression ##################################################################################################################################
# cancers with significant association between age and WGD from simple logistic regression analyses
# "OV" "SARC" "UCEC"

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("SARC")){
    model <- wgd ~ age + gender + race + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- wgd ~ age + race + figo_stage + histologic_grade
  } 
  return(model)
}

### main function
wgd_age <- function(project){
  df_tmp <- df[df$cancer_type == project,]
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  write.csv(df_tmp, paste0("Source_Data/Fig_1e_", project, ".csv", collapse = ""), row.names = FALSE)
  model <- model_selection(project)
  logit_fit <- logistf(formula = model, data=df_tmp, family = "binomial")
  summary(logit_fit)
  
  result <- broomExtra::tidy_parameters(logit_fit)
  result <- as.data.frame(result)
  p_value <- formatC(result[result$term == "age",]$p.value, format = "e", digits = 2)
  
  write.csv(result, paste0("Analysis_results/Structural_Alterations/3_Age_WGD/", project, "_multivariate_age_WGD.csv", collapse = ""), row.names = FALSE)
  
  ### Fig. 1e
  pdf(paste0("Analysis_results/Structural_Alterations/3_Age_WGD/", project, "_multivariate_age_WGD.pdf", collapse = ""), width = 3, height = 4, useDingbats=FALSE) 
  my_label <- paste0("p = ", p_value)
  p <- ggplot(df_tmp, aes(x=wgd, y=age, fill=wgd)) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(width = 0.4, fill = "white") +
    ggtitle(project) +
    xlab("WGD") +
    ylab("Age at diagnosis") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size=14,face="bold"),
          axis.title.y = element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")
  
  return(result_df)
}

cancer_types <- for_multiple_regression

df <- merge(GI, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")

results <- lapply(cancer_types, wgd_age)
results <- do.call(rbind, results)
results$q.value <- p.adjust(results$p.value, method = "BH")
results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
write.csv(results, "Analysis_results/Structural_Alterations/3_Age_WGD/Summary_multivariate_association_age_WGD.csv", row.names = FALSE)

