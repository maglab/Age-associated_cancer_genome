### PAN-Cancer Age biases in mutations (multiple logistic regression adjusting for gender, race, cancer_type)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ggplot2)
library(ggrepel)
library(reshape2)


### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to read mut burden for each cancer type and combine together
get_mut_burden <- function(project){
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id))
  mut_burden_df <- mut_burden_df[,c("Tumor_Sample_Barcode", "total")]
  mut_burden_df <- mut_burden_df[mut_burden_df$total < 1000,]
  print(paste0(project, ": ", nrow(mut_burden_df)))
  return(mut_burden_df)
}

PANCAN_mut <- lapply(projects, get_mut_burden)
PANCAN_mut <- do.call(rbind, PANCAN_mut)
dim(PANCAN_mut) # 8585 samples

###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of a gene to be mutated

# function to test association between age and mutation
test_age_mut_gene <- function(gene, mut_df, clinical){
  
  print(paste("Working on: ", gene))
  # get value for gene of interest
  values_gene <- mut_df[gene]
  
  values_gene$patient <- rownames(values_gene)
  rownames(values_gene) <- NULL
  colnames(values_gene) <- c("mut", "patient")
  
  # merge with age data
  df <- merge(values_gene, clinical, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- logistf(mut ~ age + gender + race + cancer_type , data = df, family = "binomial")
  
  summary(logit_fit)

  result <- broomExtra::tidy_parameters(logit_fit)
  
  write.csv(result, paste0("Analysis_results/Mutations/3.1_Multivariate_age_SNVs_result_each_test/PANCAN_multivariate_age_mut_", gene, ".csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- gene
  colnames(result_df) <- c("gene", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")
  result_df <- result_df[,c("gene", "estimate", "std.error", "conf.low", "conf.high", "df.error", "p.value")]

  return(result_df)
}

mut_df <- read.csv("Analysis_results/Mutations/1_table_mutations_samples/PANCAN_mutations.csv")
row.names(mut_df) <- mut_df$X
mut_df$X <- NULL
genes <- colnames(mut_df)

tmp <- do.call(rbind,lapply(genes, test_age_mut_gene, mut_df = mut_df, clinical = clinical))
rownames(tmp) <- NULL
tmp$p.value <- as.numeric(as.character(tmp$p.value))
tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
tmp$Sig <- ifelse(tmp$q.value < 0.05, TRUE, FALSE)
write.csv(tmp, "Analysis_results/Mutations/PANCAN_age_SNVs_multivariate_new.csv", row.names = FALSE)


### plot
df <- read.csv("Analysis_results/Mutations/PANCAN_age_SNVs_multivariate_new.csv")
df

### add color
my_colour <- c()
for(i in 1:nrow(df)){
  if(df$q.value[i] < 0.05 & df$estimate[i] < 0){
    my_colour <- c(my_colour, "Blue")
  } else if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
    my_colour <- c(my_colour, "Red")
  } else {
    my_colour <- c(my_colour, "Grey")
  }
}

df$my_colour <- factor(my_colour)
df$gene <- as.character(df$gene)

### Fig. 4d
pdf("Analysis_results/Mutations/PANCAN_age_SNVs_multivariate_new.pdf", width = 6, height = 4.5) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), color = my_colour, label = gene), data = df) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = c("#1d91c0", "#bdbdbd", "#a50f15")) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("PANCAN association between age and mutations") +
  xlim(c(-0.04,0.04)) +
  geom_label_repel(size = 3.5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()
