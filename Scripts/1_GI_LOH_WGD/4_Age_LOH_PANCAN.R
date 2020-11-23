# association between LOH percentage and age in PANCAN
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)

# read data
LOH <- read.table("Data/TCGA.LOH_percentage.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(LOH, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df) # n = 9678

write.csv(df, "Source_Data/Fig_1c.csv", row.names = FALSE)

### 1) Simple linear regression PANCAN ##########################################################################
lm_fit <- lm(percent_LOH ~ age, data=df)
p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$r.squared, 2)

result <- tidy(lm_fit)
write.csv(result, "Analysis_results/Structural_Alterations/4_Age_LOH/PANCAN_univariate_age_LOH_new.csv", row.names = FALSE)

### 2) Multiple linear regression PANCAN ########################################################################
# adjust on cancer type, sex, race

lm_fit <- lm(percent_LOH ~ age + cancer_type + gender + race, data=df)

p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$adj.r.squared, 2)

result <- tidy(lm_fit)
write.csv(result, "Analysis_results/Structural_Alterations/4_Age_LOH/PANCAN_multivariate_age_LOH_new.csv", row.names = FALSE)

### Fig. 1c
pdf("Analysis_results/Structural_Alterations/4_Age_LOH/PANCAN_multivariate_age_LOH_new.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
p <- ggplot(data = df, aes(x = age, y = percent_LOH)) + 
  geom_point(aes(color = cancer_type), size = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  ggtitle("Association between age and percent genome LOH") +
  xlab("Age at diagnosis") +
  ylab("percent genome LOH") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face="bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 5,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

