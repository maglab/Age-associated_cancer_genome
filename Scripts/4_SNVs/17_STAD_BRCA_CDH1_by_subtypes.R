### Subtype analysis for CDH1 mutation
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(ggplot2)
library(ggpubr)
library(rstatix)


########## STAD #################################################################################################################
# subtype information from: Liu, Y. et al. Comparative Molecular Analysis of Gastrointestinal Adenocarcinomas. Cancer Cell 33, 721-735 e728, doi:10.1016/j.ccell.2018.03.010 (2018).
project <- "STAD"

# Read mut table
mut_df <- read.csv(paste0("Analysis_results/Mutations/1_table_mutations_samples/", project, "_mutations_filter_hypermutated.csv", collapse = ""))
mut_CDH1 <- mut_df[,c("X", "CDH1")]

# Read clinical
clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))

# Read subtype data
subtypes <- read.csv("Data/TCGA_Subtypes/TCGA_GI_cancer_subtypes.csv")
subtypes <- subtypes[subtypes$TCGA.Project.Code == "STAD",]
subtypes <- subtypes[,c("TCGA.Participant.Barcode", "Molecular_Subtype")]

# merge clinical, subtype, mut
df <- merge(clin, subtypes, by.x = "patient", by.y = "TCGA.Participant.Barcode")
df <- merge(df, mut_CDH1, by.x = "patient", by.y = "X")
dim(df) # total 306 samples

df$CDH1 <- ifelse(df$CDH1 == 1, TRUE, FALSE)

df$GS <- ifelse(df$Molecular_Subtype == "GS", "GS", "Others")
df$GS <- as.factor(df$GS)

# write source data
write.csv(df, "Source_Data/Supplementary_Fig_10b.csv", row.names = FALSE)

### Plot
# subtypes and age
# Supplementary Fig 10b
# wilcoxon test
res.wilcox <- wilcox_test(age ~ GS, data = df)

write.csv(as.data.frame(res.wilcox), "Analysis_results/Mutations/7_Subtype_analysis/STAD_age_by_subtype.csv", row.names = FALSE)

pdf("Analysis_results/Mutations/7_Subtype_analysis/STAD_age_by_subtype.pdf", width = 4.5, height = 4.5) 
p <- ggviolin(df, x = "GS", y = "age", fill = "GS",
              palette = c("#edf8b1", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylab("Age at diagnosis") +
  xlab("Subtypes") +
  ggtitle("STAD: Age by subtypes") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(res.wilcox)
  )
print(p)
dev.off()

### Fisher test to investigate the enrichment of CDH1 mutations in GS subtypes
cong_table <- table(df$CDH1, df$GS)
fisher_exact_test <- fisher.test(df$CDH1, df$GS)
p_value <- fisher_exact_test$p.value
cong_table_df <- reshape2::melt(cong_table)
colnames(cong_table_df) <- c("CDH1_mutation", "Subtype", "value")

# write source data
write.csv(cong_table_df, "Source_Data/Supplementary_Fig_10a.csv", row.names = FALSE)

my_label <- paste0("p = ", signif(p_value, digits = 2))
pdf("Analysis_results/Mutations/7_Subtype_analysis/STAD_CDH1_mut_by_subtype.pdf", width = 4.5, height = 4.5)
p <- ggplot(aes(x = Subtype, y = value, fill = CDH1_mutation), data = cong_table_df) + 
  geom_bar(stat = "identity", position="fill") +
  geom_text(aes(label=value), position="fill", vjust=2) +
  ylab("proportion")  +
  labs(fill = "Mutation") +
  ggtitle("STAD: CDH1") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 6,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

########## BRCA #################################################################################################################
# subtype information from: Berger, A. C. et al. A Comprehensive Pan-Cancer Molecular Study of Gynecologic and Breast Cancers. Cancer Cell 33, 690-705 e699, doi:10.1016/j.ccell.2018.03.014 (2018).
project <- "BRCA"

# Read mut table
mut_df <- read.csv(paste0("Analysis_results/Mutations/1_table_mutations_samples/", project, "_mutations_filter_hypermutated.csv", collapse = ""))
mut_CDH1 <- mut_df[,c("X", "CDH1")]

# Read clinical
clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))

# Read subtype data
subtypes <- read.csv("Data/TCGA_Subtypes/TCGA_Gyn_and_Breast_subtypes.csv")
subtypes <- subtypes[subtypes$Tumor.Type == "BRCA",]
subtypes <- subtypes[,c("Sample.ID", "BRCA_Pathology", "BRCA_Subtype_PAM50")]

# merge clinical, subtype, mut
df <- merge(clin, subtypes, by.x = "patient", by.y = "Sample.ID")
df <- merge(df, mut_CDH1, by.x = "patient", by.y = "X")
dim(df) # total 935 samples
df$CDH1 <- ifelse(df$CDH1 == 1, TRUE, FALSE)

df_with_BRCA_Pathology <- df[!is.na(df$BRCA_Pathology),]
dim(df_with_BRCA_Pathology)   # 734 samples

df_with_BRCA_Pathology$ILC <- ifelse(df_with_BRCA_Pathology$BRCA_Pathology == "ILC", "ILC", "Others")
df_with_BRCA_Pathology$ILC <- as.factor(df_with_BRCA_Pathology$ILC)

# write source data
write.csv(df_with_BRCA_Pathology, "Source_Data/Supplementary_Fig_10d.csv", row.names = FALSE)

### Plot
# subtypes and age
# Supplementary Fig 10d
# wilcoxon test
res.wilcox <- wilcox_test(age ~ ILC, data = df_with_BRCA_Pathology)

write.csv(as.data.frame(res.wilcox), "Analysis_results/Mutations/7_Subtype_analysis/BRCA_age_by_subtype_ILC.csv", row.names = FALSE)

pdf("Analysis_results/Mutations/7_Subtype_analysis/BRCA_age_by_subtype_ILC.pdf", width = 4.5, height = 4.5) 
p <- ggviolin(df_with_BRCA_Pathology, x = "ILC", y = "age", fill = "ILC",
              palette = c("#edf8b1", "#2c7fb8"),
              add = "boxplot", add.params = list(fill = "white")) +
  ylab("Age at diagnosis") +
  xlab("Subtypes") +
  ggtitle("BRCA: Age by subtypes") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(
    subtitle = get_test_label(res.wilcox)
  )
print(p)
dev.off()

### Fisher test to investigate the enrichment of CDH1 mutations in ILC subtypes
cong_table <- table(df_with_BRCA_Pathology$CDH1, df_with_BRCA_Pathology$ILC)
fisher_exact_test <- fisher.test(df_with_BRCA_Pathology$CDH1, df_with_BRCA_Pathology$ILC)
p_value <- fisher_exact_test$p.value
cong_table_df <- reshape2::melt(cong_table)
colnames(cong_table_df) <- c("gene_mutation", "Subtype", "value")

# write source data
write.csv(cong_table_df, "Source_Data/Supplementary_Fig_10c.csv", row.names = FALSE)

my_label <- paste0("p = ", signif(p_value, digits = 2))

pdf("Analysis_results/Mutations/7_Subtype_analysis/BRCA_CDH1_mut_by_subtype_ILC.pdf", width = 4.5, height = 4.5)
p <- ggplot(aes(x = Subtype, y = value, fill = gene_mutation), data = cong_table_df) + 
  geom_bar(stat = "identity", position="fill") +
  geom_text(aes(label=value), position="fill", vjust=2) +
  ylab("proportion")  +
  labs(fill = "Mutation") +
  ggtitle("BRCA: CDH1") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 6,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

