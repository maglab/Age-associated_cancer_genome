### association between age and overall SCNA score in LUAD by smoking status
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(ggpubr)
library(rstatix)

### Read data
# read SCNA score
SCNA <- read.csv("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/LUAD_normalised_SCNA_score.csv")
colnames(SCNA) <- c("patient", "rank_norm_focal", "rank_norm_chrom", "rank_norm_arm", "sum_score")
dim(SCNA)

# read clinical
clin <- read.csv("Data/clinical_XML_interest/TCGA-LUAD_clinical_XML_interest.csv")

# merge
df <- merge(SCNA, clin, by.x = "patient", by.y = "patient")
head(df)
df$smoking_history <- ordered(df$smoking_history, levels = c("Current smoker", "Current reformed smoker", "Lifelong Non-smoker"))
df_no_NA <- df[!is.na(df$smoking_history),]

### linear regression SCNA and age (only Lifelong Non-smoker)
df_no_smokers <- df_no_NA[df_no_NA$smoking_history == "Lifelong Non-smoker",]

# write source data
write.csv(df_no_smokers, "Source_Data/Supplementary_Fig_4a.csv", row.names = FALSE)

lm_fit <- lm(sum_score ~ age + gender + race + pathologic_stage, data=df_no_smokers) 
summary(lm_fit)

p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$adj.r.squared, 4)

my_label <- paste0("adj.R-squared = ", r_squared, "\np = ", p_value)
pdf("Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Lifelong_Non-smoker.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(data = df_no_smokers, aes(x = age, y = sum_score)) + 
  geom_point(color='black') +
  geom_smooth(method = "lm") +
  ggtitle("LUAD (Only Lifelong Non-Smoker)") +
  xlab("Age at diagnosis") +
  ylab("Overall SCNA score") +
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
write.csv(result, "Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Lifelong_Non-smoker.csv", row.names = FALSE)

### linear regression SCNA and age (only Current smokers)
df_smokers <- df_no_NA[df_no_NA$smoking_history == "Current smoker",]

# write source data
write.csv(df_smokers, "Source_Data/Supplementary_Fig_4b.csv", row.names = FALSE)

lm_fit <- lm(sum_score ~ age + gender + race + pathologic_stage, data=df_smokers) 
summary(lm_fit)

p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$adj.r.squared, 4)

my_label <- paste0("adj.R-squared = ", r_squared, "\np = ", p_value)
pdf("Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Current_smoker.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(data = df_smokers, aes(x = age, y = sum_score)) + 
  geom_point(color='black') +
  geom_smooth(method = "lm") +
  #ggtitle(paste0("Association between age and GI score in ", project, collapse = "")) +
  ggtitle("LUAD (Only Current Smoker)") +
  xlab("Age at diagnosis") +
  ylab("Overall SCNA score") +
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

write.csv(result, "Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Current_smoker.csv", row.names = FALSE)

### linear regression SCNA and age (only Current reformed smokers)
df_reformed_smokers <- df_no_NA[df_no_NA$smoking_history == "Current reformed smoker",]

# write source data
write.csv(df_reformed_smokers, "Source_Data/Supplementary_Fig_4c.csv", row.names = FALSE)

lm_fit <- lm(sum_score ~ age + gender + race + pathologic_stage, data=df_reformed_smokers) 
summary(lm_fit)

p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$adj.r.squared, 4)

my_label <- paste0("adj.R-squared = ", r_squared, "\np = ", p_value)
pdf("Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Current_reformed_smoker.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(data = df_reformed_smokers, aes(x = age, y = sum_score)) + 
  geom_point(color='black') +
  geom_smooth(method = "lm") +
  #ggtitle(paste0("Association between age and GI score in ", project, collapse = "")) +
  ggtitle("LUAD (Only Current Reformed Smoker)") +
  xlab("Age at diagnosis") +
  ylab("Overall SCNA score") +
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

write.csv(result, "Analysis_results/CNAs/2_SCNA_scores/LUAD_multivariate_age_overall_SCNA_Current_reformed_smoker.csv", row.names = FALSE)

