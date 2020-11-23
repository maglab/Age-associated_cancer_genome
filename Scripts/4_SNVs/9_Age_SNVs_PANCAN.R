### PAN-Cancer Age biases in mutations (multiple logistic regression adjusting for gender, race, cancer_type)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)

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
  colnames(result_df) <- c("gene", "estimate", "std.error", "conf.low", "conf.high", "statistic", "df.error", "p.value")
  result_df <- result_df[,c("gene", "estimate", "std.error", "conf.low", "conf.high", "statistic", "df.error", "p.value")]
  
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

# write source data
write.csv(df, "Source_Data/Fig_4e.csv", row.names = FALSE)

### Fig. 4e
# plot
pdf("Analysis_results/Mutations/PANCAN_age_SNVs_multivariate_new.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), color = my_colour, label = gene), data = df) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("#1d91c0", "#bdbdbd", "#a50f15")) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("PANCAN: Age and mutations") +
  xlim(c(-0.041,0.041)) +
  geom_label_repel(size = 4.5) +
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
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()
