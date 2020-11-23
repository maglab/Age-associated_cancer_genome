### Age distribution of LUAD, ESCA, LIHC
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(ggpubr)
library(rstatix)

########## LUAD and smoking #####################################################################################################
### Read data
samples_of_interest <- read.csv("Data/samples_in_ASCAT_and_seg.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# read LOH
LOH <- read.table("Data/TCGA.LOH_not_exclude_aneuploidy.clean.txt", header = TRUE)

# merge data
df <- merge(LOH, purity[, c("patient", "purity")], by.x = "patient", by.y = "patient")
df <- df[df$patient %in% samples_of_interest$patient,]
dim(df)

# read clinical
clin <- read.csv("Data/clinical_XML_interest/TCGA-LUAD_clinical_XML_interest.csv")

# merge
df_LUAD <- merge(df, clin, by.x = "patient", by.y = "patient")
head(df_LUAD)
df_LUAD$smoking_history <- ordered(df_LUAD$smoking_history, levels = c("Current smoker", "Current reformed smoker", "Lifelong Non-smoker"))

### Plot
# smoking status and age
# kruskal test
df_no_NA <- df_LUAD[!is.na(df_LUAD$smoking_history),]
res.kruskal <- kruskal_test(age ~ smoking_history, data = df_no_NA)
pwc <- dunn_test(age ~ smoking_history, data = df_no_NA, p.adjust.method = "bonferroni") # pairwise comparisons using dunn test
pwc$p.adj.format <- formatC(pwc$p.adj, format = "e", digits = 2)  # format p-values to scientific notation

write.csv(as.data.frame(pwc), "Analysis_results/Structural_Alterations/4_Age_LOH/LUAD_age_by_smoking.csv", row.names = FALSE)

# write source data
write.csv(df_no_NA, "Source_Data/Supplementary_Fig_2a.csv", row.names = FALSE)

# plot
pdf("Analysis_results/Structural_Alterations/4_Age_LOH/LUAD_age_by_smoking.pdf", width = 4.5, height = 6, useDingbats=FALSE)
p <- ggviolin(df_no_NA, x = "smoking_history", y = "age", fill = "smoking_history",
              palette = c("#edf8b1", "#7fcdbb", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylim(20, 110) +
  ylab("Age at diagnosis") +
  xlab("Smoking status") +
  ggtitle("LUAD: Age by smoking status") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  stat_pvalue_manual(pwc[pwc$group1 == "Current smoker" | pwc$group1 == "Current reformed smoker",], label = "p.adj.format",
                     hide.ns = FALSE, y.position = c(max(df_no_NA$age) + 5, max(df_no_NA$age) + 15, max(df_no_NA$age) + 10)) +
  labs(
    subtitle = get_test_label(res.kruskal)
  )
print(p)
dev.off()




########## ESCA and smoking #####################################################################################################
# read clinical
clin <- read.csv("Data/clinical_XML_interest/TCGA-ESCA_clinical_XML_interest.csv")

# merge
df_ESCA <- merge(df, clin, by.x = "patient", by.y = "patient")
head(df_ESCA)
dim(df_ESCA)

### smoking
df_ESCA$smoking_history <- ordered(df_ESCA$smoking_history, levels = c("Current smoker", "Current reformed smoker", "Lifelong Non-smoker"))

# smoking status and age
# kruskal test
df_no_NA <- df_ESCA[!is.na(df_ESCA$smoking_history),]
res.kruskal <- kruskal_test(age ~ smoking_history, data = df_no_NA)
pwc <- dunn_test(age ~ smoking_history, data = df_no_NA, p.adjust.method = "bonferroni") # pairwise comparisons using dunn test
pwc$p.adj.format <- formatC(pwc$p.adj, format = "e", digits = 2)  # format p-values to scientific notation

write.csv(as.data.frame(pwc), "Analysis_results/Structural_Alterations/4_Age_LOH/ESCA_age_by_smoking.csv", row.names = FALSE)

# write source data
write.csv(df_no_NA, "Source_Data/Supplementary_Fig_2b.csv", row.names = FALSE)

pdf("Analysis_results/Structural_Alterations/4_Age_LOH/ESCA_age_by_smoking.pdf", width = 4.5, height = 6, useDingbats=FALSE) 
p <- ggviolin(df_no_NA, x = "smoking_history", y = "age", fill = "smoking_history",
              palette = c("#edf8b1", "#7fcdbb", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylab("Age at diagnosis") +
  xlab("Smoking status") +
  ggtitle("ESCA: Age by smoking status") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  stat_pvalue_manual(pwc[pwc$group1 == "Current smoker" | pwc$group1 == "Current reformed smoker",], label = "p.adj.format",
                     hide.ns = FALSE, y.position = c(max(df_no_NA$age) + 5, max(df_no_NA$age) + 15, max(df_no_NA$age) + 10)) +
  labs(
    subtitle = get_test_label(res.kruskal)
  )
print(p)
dev.off()

########## ESCA and race ########################################################################################################
# race and age
# kruskal test
df_no_NA <- df_ESCA[!is.na(df_ESCA$race),]
res.kruskal <- kruskal_test(age ~ race, data = df_no_NA)
pwc <- dunn_test(age ~ race, data = df_no_NA, p.adjust.method = "bonferroni") # pairwise comparisons using dunn test
pwc$p.adj.format <- formatC(pwc$p.adj, format = "e", digits = 2)  # format p-values to scientific notation

write.csv(as.data.frame(pwc), "Analysis_results/Structural_Alterations/4_Age_LOH/ESCA_age_by_race.csv", row.names = FALSE)

# write source data
write.csv(df_no_NA, "Source_Data/Supplementary_Fig_2c.csv", row.names = FALSE)

pdf("Analysis_results/Structural_Alterations/4_Age_LOH/ESCA_age_by_race.pdf", width = 4.5, height = 6, useDingbats=FALSE)
p <- ggviolin(df_no_NA, x = "race", y = "age", fill = "race",
              palette = c("#edf8b1", "#7fcdbb", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylab("Age at diagnosis") +
  xlab("Race") +
  ggtitle("ESCA: Age by race") +
  scale_x_discrete(labels=c("Asian", "Black/African American", "White")) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  stat_pvalue_manual(pwc[pwc$group1 == "ASIAN" | pwc$group1 == "BLACK OR AFRICAN AMERICAN",], label = "p.adj.format",
                     hide.ns = FALSE, y.position = c(max(df_no_NA$age) + 5, max(df_no_NA$age) + 15, max(df_no_NA$age) + 10)) +
  labs(
    subtitle = get_test_label(res.kruskal)
  )
print(p)
dev.off()

########## LIHC and histologic_grade ############################################################################################
# read clinical
clin <- read.csv("Data/clinical_XML_interest/TCGA-LIHC_clinical_XML_interest.csv")

# merge
df_LIHC <- merge(df, clin, by.x = "patient", by.y = "patient")
head(df_LIHC)
dim(df_LIHC)

### Plot
# histologic grade and age
# kruskal test
df_no_NA <- df_LIHC[!is.na(df_LIHC$histologic_grade),]
res.kruskal <- kruskal_test(age ~ histologic_grade, data = df_no_NA)
pwc <- dunn_test(age ~ histologic_grade, data = df_no_NA, p.adjust.method = "bonferroni") # pairwise comparisons using dunn test
pwc$p.adj.format <- formatC(pwc$p.adj, format = "e", digits = 2)  # format p-values to scientific notation

write.csv(as.data.frame(pwc), "Analysis_results/Structural_Alterations/4_Age_LOH/LIHC_age_by_histologic_grade.csv", row.names = FALSE)

# write source data
write.csv(df_no_NA, "Source_Data/Supplementary_Fig_2d.csv", row.names = FALSE)

pdf("Analysis_results/Structural_Alterations/4_Age_LOH/LIHC_age_by_histologic_grade.pdf", width = 4.5, height = 6, useDingbats=FALSE) 
p <- ggviolin(df_no_NA, x = "histologic_grade", y = "age", fill = "histologic_grade",
              palette = c("#edf8b1", "#7fcdbb", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylab("Age at diagnosis") +
  xlab("Histologic grade") +
  ggtitle("LIHC: Age by histologic grade") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  stat_pvalue_manual(pwc[pwc$group1 == "G1" | pwc$group1 == "G2",], label = "p.adj.format",
                     hide.ns = FALSE, y.position = c(max(df_no_NA$age) + 5, max(df_no_NA$age) + 15, max(df_no_NA$age) + 10)) +
  labs(
    subtitle = get_test_label(res.kruskal)
  )
print(p)
dev.off()


