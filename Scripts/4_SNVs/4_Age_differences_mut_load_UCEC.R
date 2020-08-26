### UCEC Age biases in mutations

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(reshape2)

### read data
project <- "UCEC"
clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# mut burden
mut_burden <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
mut_burden$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden$Tumor_Sample_Barcode), clean_id))

mut_burden <- mut_burden[,c("Tumor_Sample_Barcode", "total")]
mut_burden$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden$Tumor_Sample_Barcode), clean_id))

df <- merge(clin, mut_burden, by.x = "patient", by.y = "Tumor_Sample_Barcode")
df$mut_burden_group <- ifelse(df$total > 1000, ">1000", "<=1000")
df$age_group <- ifelse(df$age > 50, ">50", "<=50")
cong_table <- table(df$mut_burden_group, df$age_group)
fisher_exact_test <- fisher.test(df$mut_burden_group, df$age_group)
p_value <- round(fisher_exact_test$p.value, 4)
cong_table_df <- reshape2::melt(cong_table)
colnames(cong_table_df) <- c("Mutational_burden", "Age", "value")

my_label <- paste0("p = ", p_value)

### Fig. 4a
pdf("Analysis_results/Mutations/5_UCEC_low_burden_age_SNVs/Age_and_mut_burden_stacked.pdf", width = 4, height = 4)
p <- ggplot(aes(x = Age, y = value, fill = Mutational_burden), data = cong_table_df) + 
  geom_bar(stat = "identity", position="fill") +
  geom_text(aes(label=value), position="fill", vjust=2) +
  ylab("proportion")  +
  ggtitle("UCEC mutations and age") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        panel.background = element_blank(),
        #panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 5,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

