# association between WGD status and age in PANCAN
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)
library(logistf)

# read data
samples_of_interest <- read.csv("Data/samples_in_ASCAT_and_seg.csv")

GI <- read.table("Data/TCGA.giScores.wgd.txt", header = TRUE)
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# merge data
df <- merge(GI, clinical[,c("patient", "age", "gender", "race")], by.x = "patient", by.y = "patient")
df <- merge(df, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
dim(df) # 9678

### 1) Simple logistic regression PANCAN ##########################################################################
logit_fit <- glm(wgd ~ age , data = df, family = "binomial")
summary(logit_fit)
result <- tidy(logit_fit)
p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
result <- as.data.frame(parameters::parameters(logit_fit))
write.csv(result, "Analysis_results/Structural_Alterations/3_Age_WGD/PANCAN_univariate_age_WGD.csv", row.names = FALSE)

pdf("Analysis_results/Structural_Alterations/3_Age_WGD/PANCAN_univariate_age_WGD.pdf", width = 3, height = 4) 
my_label <- paste0("p = ", p_value)
p <- ggplot(df, aes(x=wgd, y=age, fill=wgd)) + 
  geom_violin(trim = FALSE, scale = "width") + 
  geom_boxplot(width = 0.4, fill = "white") +
  #geom_jitter(aes(color = cancer_type)) +
  #ggtitle("Association between age and whole genome duplication") +
  ggtitle("PANCAN") +
  xlab("Whole Genome Duplication") +
  ylab("Age at diagnosis") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, 
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

### 2) Multiple logistic regression PANCAN ########################################################################
# adjust on cancer type, sex, race
logit_fit <- glm(wgd ~ age + cancer_type + gender + race, data = df, family = "binomial")
summary(logit_fit)

result <- tidy(logit_fit)
p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
result <- as.data.frame(parameters::parameters(logit_fit))

write.csv(result, "Analysis_results/Structural_Alterations/3_Age_WGD/PANCAN_multivariate_age_WGD.csv", row.names = FALSE)

pdf("Analysis_results/Structural_Alterations/3_Age_WGD/PANCAN_multivariate_age_WGD.pdf", width = 3, height = 4) 
my_label <- paste0("p = ", p_value)
p <- ggplot(df, aes(x=wgd, y=age, fill=wgd)) +
  geom_violin(trim = FALSE, scale = "width") + 
  geom_boxplot(width = 0.4, fill = "white") +
  ggtitle("PANCAN") +
  xlab("Whole Genome Duplication") +
  ylab("Age at diagnosis") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, 
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

