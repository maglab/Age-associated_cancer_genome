# association between percent genomic LOH and age in each cancer type
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)

# read data
LOH <- read.table("Data/TCGA.LOH_not_exclude_aneuploidy.clean.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(LOH, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df) # n = 9678

cancer_types <- sort(unique(as.character(df$cancer_type)))

# create a function to test the association between percent genomic LOH and age
LOH_age <- function(cancer_type, df){
  df <- df[df$cancer_type == cancer_type,]
  print(paste0(cancer_type, " - num sample: ", nrow(df), collapse = ""))
  
  median_LOH <- median(df$percent_LOH)
  
  lm_fit <- lm(percent_LOH ~ age, data=df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared, 2)
  
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  result <- tidy(lm_fit)
  write.csv(result, paste0("Analysis_results/Structural_Alterations/4_Age_LOH/", cancer_type, "_univariate_age_LOH_new.csv", collapse = ""), row.names = FALSE)
  
  pdf(paste0("Analysis_results/Structural_Alterations/4_Age_LOH/", cancer_type, "_univariate_age_LOH_new.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE) 
  my_label <- paste0("R-squared = ", r_squared, "\np = ", p_value)
  p <- ggplot(data = df, aes(x = age, y = percent_LOH)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    ggtitle(cancer_type) +
    xlab("Age at diagnosis") +
    ylab("percent genomic LOH") +
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
  result_df$median_LOH <- median_LOH
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value", "median_LOH")
  result_df <- result_df[,c("cancer_type", "estimate", "std.error", "statistic", "median_LOH", "p.value")]
  return(result_df)
}

LOH_results <- lapply(cancer_types, LOH_age, df = df)
LOH_results <- do.call(rbind, LOH_results)

LOH_results$estimate <- as.numeric(as.character(LOH_results$estimate))
LOH_results$median_LOH <- as.numeric(as.character(LOH_results$median_LOH))
LOH_results$p.value <- as.numeric(as.character(LOH_results$p.value))

LOH_results$q.value <- p.adjust(LOH_results$p.value, method = "BH")
LOH_results$Sig <- ifelse(LOH_results$q.value < 0.05, TRUE, FALSE)

write.csv(LOH_results, "Analysis_results/Structural_Alterations/4_Age_LOH/Summary_univariate_association_age_LOH_new.csv", row.names = FALSE)

for_multiple_regression <- LOH_results[LOH_results$Sig == TRUE,]$cancer_type

### Multiple linear regression ##################################################################################################################################
# cancers with significant association between age and percent genomic LOH from simple linear regression analyses
# "CESC" "ESCA" "LGG"  "LIHC" "LUAD" "PRAD" "SARC" "THCA" "UCEC" "UVM"
# each cancer type will be analysed separately, using different model

model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- percent_LOH ~ age + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- percent_LOH ~ age + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- percent_LOH ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- percent_LOH ~ age + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- percent_LOH ~ age + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- percent_LOH ~ age + gender + race
  } else if(project %in% c("HNSC")){
    model <- percent_LOH ~ age + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- percent_LOH ~ age + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- percent_LOH ~ age + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- percent_LOH ~ age + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- percent_LOH ~ age + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- percent_LOH ~ age + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- percent_LOH ~ age + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- percent_LOH ~ age + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- percent_LOH ~ age + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- percent_LOH ~ age + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- percent_LOH ~ age + gender + pathologic_stage
  }
  return(model)
}

### main function
LOH_age <- function(project){
  df_tmp <- df[df$cancer_type == project,]
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  df_tmp <- merge(df_tmp, clin, by.x = "patient", by.y = "patient")
  head(df_tmp)
  
  write.csv(df_tmp, paste0("Source_Data/Supplementary_Fig_1b_", project, ".csv"), row.names = FALSE)
  
  model <- model_selection(project)
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  r_squared <- round(summary(lm_fit)$adj.r.squared,2)
  
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  pdf(paste0("Analysis_results/Structural_Alterations/4_Age_LOH/", project, "_multivariate_age_LOH_new.pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE) 
  p <- ggplot(data = df_tmp, aes(x = age, y = percent_LOH)) + 
    geom_point(color='black') +
    geom_smooth(method = "lm") +
    ggtitle(project) +
    xlab("Age at diagnosis") +
    ylab("percent genomic LOH") +
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
  write.csv(result, paste0("Analysis_results/Structural_Alterations/4_Age_LOH/", project, "_multivariate_age_LOH_new.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[2,])
  result_df$term <- project
  colnames(result_df) <- c("cancer_type", "estimate", "std.error", "statistic", "p.value")
  return(result_df)
}

cancer_types <- for_multiple_regression

df <- merge(LOH, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")

result <- lapply(cancer_types, LOH_age)
result <- do.call(rbind, result)

# adjust p-value
result$q.value <- p.adjust(result$p.value, method = "BH")
result$Sig <- ifelse(result$q.value < 0.05, TRUE, FALSE)

write.csv(result, "Analysis_results/Structural_Alterations/4_Age_LOH/Summary_multivariate_association_age_LOH_new.csv", row.names = FALSE)
write.csv(result, "Source_Data/Fig_1d.csv", row.names = FALSE)

### Fig. 1d Volcano plot ##################################################################################
### add information on whether each cancer still showed a significant association after multiple regression
LOH_results$LOH_multivariate <- ifelse(LOH_results$cancer_type %in% result[result$Sig == TRUE,]$cancer_type, TRUE, FALSE)

### add color
my_colour <- c()
for(i in 1:nrow(LOH_results)){
  if(LOH_results$LOH_multivariate[i] == TRUE & LOH_results$estimate[i] < 0){
    my_colour <- c(my_colour, "Blue")
  } else if(LOH_results$LOH_multivariate[i] == TRUE & LOH_results$estimate[i] > 0){
    my_colour <- c(my_colour, "Red")
  } else if(LOH_results$Sig[i] == TRUE & LOH_results$LOH_multivariate[i] == FALSE){
    my_colour <- c(my_colour, "Black")
  }else {
    my_colour <- c(my_colour, "Grey")
  }
}

LOH_results$my_colour <- factor(my_colour)

# plot
pdf("Analysis_results/Structural_Alterations/4_Age_LOH/Summary_multivariate_association_age_LOH_new.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), size = median_LOH, 
                color = my_colour, label = cancer_type), data = LOH_results) +
  geom_point() +
  scale_color_manual(values = c("#000000", "#1d91c0", "#bdbdbd", "#a50f15")) +
  scale_size_continuous(range = c(1, 4)) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Age and percent genome LOH by cancer type") +
  xlim(c(-0.35,0.35)) +
  geom_label_repel(size = 4) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
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
# - Blue: significant multiple linear regression and negative coefficient
# Size corresponses to median percent LOH


