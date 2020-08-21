# association between GI score and age in PANCAN
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)

# read data
GI <- read.table("Data/TCGA.giScores.wgd.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(GI, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df) # n = 9678

### Multiple linear regression PANCAN ########################################################################
# adjust on cancer type, sex, race

lm_fit <- lm(gi ~ age + cancer_type + gender + race, data=df)

p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
r_squared <- round(summary(lm_fit)$adj.r.squared, 2)

result <- tidy(lm_fit)
write.csv(result, "Analysis_results/Structural_Alterations/1_Age_GIscore/PANCAN_multivariate_age_GIscore.csv", row.names = FALSE)

### Fig 1a
pdf("Analysis_results/Structural_Alterations/1_Age_GIscore/PANCAN_multivariate_age_GIscore.pdf", width = 8, height = 6) 
my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
p <- ggplot(data = df, aes(x = age, y = gi)) + 
  geom_point(aes(color = cancer_type), size = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  ggtitle("Association between age and GI score") +
  xlab("Age at diagnosis") +
  ylab("GI score") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face="bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 5,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()
