# association between GI score and age in each cancer type
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)
library(ggrepel)

# read data
GI <- read.table("Data/TCGA.giScores.wgd.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(GI, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df)

cancer_types <- sort(as.character(unique(df$cancer_type)))

# create a function to test the association between GI score and age
GI_age <- function(cancer_type){
  df <- df[df$cancer_type == cancer_type,]
  print(paste0(cancer_type, " - num sample: ", nrow(df), collapse = ""))
  
  median_GI <- median(df$gi)
  
  lm_fit <- lm(gi ~ age, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared, 2)
  
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  result <- tidy(lm_fit)
  write.csv(result, paste0("Analysis_results/Structural_Alterations/1_Age_GIscore/", cancer_type, "_univariate_age_GIscore.csv", collapse = ""), row.names = FALSE)
  
  my_label <- paste0("R-squared = ", r_squared, "\np = ", p_value)
  pdf(paste0("Analysis_results/Structural_Alterations/1_Age_GIscore/", cancer_type, "_univariate_age_GIscore.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE) 
  p <- ggplot(data = df, aes(x = age, y = gi)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    ggtitle(cancer_type) +
    xlab("Age at diagnosis") +
    ylab("GI score") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 6,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- cancer_type
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  result_df$median_GI <- median_GI
  return(result_df)
}

GI_result <- lapply(cancer_types, GI_age)
GI_result <- do.call(rbind, GI_result)

GI_result$estimate <- as.numeric(as.character(GI_result$estimate))
GI_result$median_GI <- as.numeric(as.character(GI_result$median_GI))
GI_result$p.value <- as.numeric(as.character(GI_result$p.value))

GI_result$q.value <- p.adjust(GI_result$p.value, method = "BH")
GI_result$Sig <- ifelse(GI_result$q.value < 0.05, TRUE, FALSE)
GI_result <- GI_result[,c("cancer_type", "median_GI", "estimate", "std.error", "statistic", "p.value", "q.value", "Sig")]

write.csv(GI_result, "Analysis_results/Structural_Alterations/1_Age_GIscore/Summary_Association_Age_univariate_GIscore.csv", row.names = FALSE)

for_multiple_regression <- GI_result[GI_result$Sig == TRUE,]$cancer_type

### Multiple linear regression ##################################################################################################################################
# cancers with significant association between age and GI score from simple linear regression analysis
# BLCA, CESC, LGG, LUAD, OV, PRAD, SARC, TGCT, THCA, UCEC
# each cancer type will be analysed separately, using different model

# sub-function to select a regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- gi ~ age + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- gi ~ age + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- gi ~ age + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- gi ~ age + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- gi ~ age + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- gi ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- gi ~ age + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- gi ~ age + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- gi ~ age + gender + race
  } else if(project %in% c("HNSC")){
    model <- gi ~ age + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- gi ~ age + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- gi ~ age + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- gi ~ age + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- gi ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- gi ~ age + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- gi ~ age + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- gi ~ age + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- gi ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- gi ~ age + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- gi ~ age + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- gi ~ age + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- gi ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- gi ~ age + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- gi ~ age + gender + pathologic_stage
  }
  return(model)
}

### main function
gi_age <- function(project){
  df_tmp <- df[df$cancer_type == project,]
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  
  write.csv(df_tmp, paste0("Source_Data/Supplementary_Fig_1a_", project, ".csv"), row.names = FALSE)
  
  model <- model_selection(project)
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$adj.r.squared,2)
  
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)

  ### Supplementary Fig. 1a
  pdf(paste0("Analysis_results/Structural_Alterations/1_Age_GIscore/", project, "_multivariate_age_GIscore.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE) 
  p <- ggplot(data = df_tmp, aes(x = age, y = gi)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("GI score") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 6,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  
  result <- tidy(lm_fit)
  write.csv(result, paste0("Analysis_results/Structural_Alterations/1_Age_GIscore/", project, "_multivariate_age_GIscore.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  return(result_df)
}

cancer_types <- for_multiple_regression

df <- merge(GI, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")

results <- lapply(cancer_types, gi_age)
results <- do.call(rbind, results)

# adjust p-value
results$q.value <- p.adjust(results$p.value, method = "BH")

results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
write.csv(results, "Analysis_results/Structural_Alterations/1_Age_GIscore/Summary_Association_Age_multivariate_GIscore.csv", row.names = FALSE)
write.csv(results, "Source_Data/Fig_1b.csv", row.names = FALSE)

### Fig. 1b Volcano plot ##################################################################################
### add information on whether each cancer still showed a significant association after multiple regression
GI_result$GI_multivariate <- ifelse(GI_result$cancer_type %in% results[results$Sig == TRUE,]$cancer_type, TRUE, FALSE)
GI_result$cancer_type <- as.character(GI_result$cancer_type)

### add color
my_colour <- c()
for(i in 1:nrow(GI_result)){
  if(GI_result$q.value[i] < 0.05 & GI_result$GI_multivariate[i] == FALSE){
    my_colour <- c(my_colour, "Black")
  } else if(GI_result$GI_multivariate[i] == TRUE & GI_result$estimate[i] > 0){
    my_colour <- c(my_colour, "Red")
  } else {
    my_colour <- c(my_colour, "Grey")
  }
}

GI_result$my_colour <- factor(my_colour)

pdf("Analysis_results/Structural_Alterations/1_Age_GIscore/Summary_multivariate_age_GIscore.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), size = median_GI, 
                color = my_colour, label = cancer_type), data = GI_result) +
  geom_point() +
  scale_color_manual(values = c("#000000", "#bdbdbd", "#a50f15")) +
  scale_size_continuous(range = c(1, 4)) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Age and GI score by cancer type") +
  xlim(c(-0.0085,0.0085)) +
  geom_label_repel(size = 4) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5)+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()


# here in the volcano plot are the plot between 
# - X axis: simple linear regression coefficient
# - Y axis: simple linear regression -log10(adjusted p value by BH method)
# Colour represents significant
# - Grey: not significant
# - Red: significant in multiple linear regression and positive coefficient
# - Black: significant in simple but not in multiple regression analyses
# Size corresponses to median GI score.

