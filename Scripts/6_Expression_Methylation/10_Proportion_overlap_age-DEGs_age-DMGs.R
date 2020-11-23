### barchart between age-DEGs and age-DMGs

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)

### projects with > 150 age-DEGs and > 150 age-DMGs
projects <- c("LGG", "BRCA", "UCEC", "ESCA", "KIRP", "OV", "LIHC", "LAML", "SKCM", "PRAD")

### function to plot proportion
get_num_genes <- function(project){
  
  # read files
  methylation <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  met_up <- as.character(methylation[methylation$Sig == TRUE & methylation$estimate > 0,]$gene)
  met_down <- as.character(methylation[methylation$Sig == TRUE & methylation$estimate < 0,]$gene)
  DMGs <- c(met_up, met_down)
  
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  exp_up <- as.character(expression[expression$Sig == TRUE & expression$estimate > 0,]$gene)
  exp_down <- as.character(expression[expression$Sig == TRUE & expression$estimate < 0,]$gene)
  DEGs <- c(exp_up, exp_down)
  
  overlap_genes <- DMGs[DMGs %in% DEGs]
  
  met_overlap_genes <- methylation[methylation$gene %in% overlap_genes,]
  met_overlap_genes$direction_methy <- ifelse(met_overlap_genes$estimate > 0, "Methy_up", "Methy_down")
  met_overlap_genes <- met_overlap_genes[,c("gene", "direction_methy")]
  
  exp_overlap_genes <- expression[expression$gene %in% overlap_genes,]
  exp_overlap_genes$direction_expr <- ifelse(exp_overlap_genes$estimate > 0, "Expr_up", "Expr_down")
  exp_overlap_genes <- exp_overlap_genes[,c("gene", "direction_expr")]
  
  df <- merge(met_overlap_genes, exp_overlap_genes, by.x = "gene", by.y = "gene")
  df$condition <- paste(df$direction_methy, "; ", df$direction_expr, sep = "")
  
  df <- as.data.frame(table(df$condition))
  df$project <- rep(project, nrow(df))
  colnames(df) <- c("Condition", "Freq", "Project")
  return(df)
}

my_df <- lapply(projects, get_num_genes)

my_df <- do.call(rbind, my_df)
my_df <- my_df[,c("Project", "Condition", "Freq")]

write.csv(my_df, "Analysis_results/Methylation/Overlap_age-DMGs_age-DEGs.csv", row.names = FALSE)

# write source data
write.csv(my_df, "Source_Data/Fig_6c.csv", row.names = FALSE)

### Fig. 6c
pdf("Analysis_results/Methylation/Overlap_age-DMGs_age-DEGs.pdf", width = 8, height = 4) 
p <- ggplot(data = my_df, aes(x = Project, y = Freq, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("proportion")  +
  ggtitle("Overlap genes between age-DMGs and age-DEGs") +
  scale_fill_manual(values = c("#0571b0", "#f4a582", "#92c5de", "#ca0020")) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_text(size = 12,face="bold"),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()

