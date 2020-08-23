### Summarise genes with sig association between SCNAs and gene expression
### including cancer driver genes
### Fig. 3d

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(reshape2)
library(ggplot2)

files <- list.files(path = "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/SCNAs_gene_expression", 
                    pattern = "*.csv",
                    full.names = TRUE)

# remove LIHC gain because it does not contain any gene
files <- files[files != "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/SCNAs_gene_expression/LIHC_CNAs_Gene_Expression_gain.csv"]

# function to keep only sig genes and add cancer_type, GainOrLoss
get_sig_genes <- function(file){
  
  print(paste0("Working on: ", file))
  
  short_file_name <- unlist(strsplit(file, split = "/"))[6]
  
  # check gain or del
  if(grepl(x = short_file_name, pattern = "gain")){
    GainOrLoss <- "gain"
  } else if (grepl(x = short_file_name, pattern = "loss")){
    GainOrLoss <- "loss"
  }
  
  tmp <- read.csv(file)
  cancer_type <- unlist(strsplit(short_file_name, split = "_"))[1]
  tmp$cancer_type <- rep(cancer_type, nrow(tmp))
  tmp$GainOrLoss <- rep(GainOrLoss, nrow(tmp))
  colnames(tmp) <- c("gene", colnames(tmp)[2:ncol(tmp)])
  tmp <- tmp[,c("cancer_type", "gene", "estimate", "statistic", "p.value", "parameter", 
                "conf.low", "conf.high", "method", "alternative", "adj_p", "sig", "GainOrLoss")]
  #tmp <- tmp[tmp$sig == TRUE,]
  tmp <- tmp[complete.cases(tmp),]
  return(tmp)
}

df <- lapply(files, get_sig_genes)
df <- do.call(rbind, df)

### cancer driver genes
cancer_genes <- as.character(unlist(read.table("Data/cancer_genes.txt")))

df$Cancer_driver <- ifelse(df$gene %in% cancer_genes, TRUE, FALSE)

write.csv(df, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_corr_CNAs_gene_expression.csv", row.names = FALSE)

# only sig
df <- df[df$sig == TRUE,]
write.csv(df, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_sig_corr_CNAs_gene_expression.csv", row.names = FALSE)

### Plot
### read data
gene_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_sig_corr_CNAs_gene_expression.csv")
gain <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_gain.csv")
loss <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_del.csv")

### read cancer-related genes

gene_df <- gene_df[gene_df$Cancer_driver == TRUE,]

### get direction (increase or decrease with age)
get_direction <- function(i, gene_df, gain, loss){
  
  tmp <- gene_df[i,]
  
  cancer_type <- as.character(tmp$cancer_type)
  gene <- as.character(tmp$gene)
  print(paste0("Working on: ", cancer_type, ", ", gene))
  GainOrLoss <- as.character(tmp$GainOrLoss)
  
  if (GainOrLoss == "gain"){
    my_df <- gain
  } else if (GainOrLoss == "loss"){
    my_df <- loss
  }
  
  my_df <- my_df[my_df$cancer_type == cancer_type, ]
  my_df <- my_df[grepl(gene, my_df$genes),]
  
  if(nrow(my_df) > 1){
    my_df <- my_df[my_df$q.value == min(my_df$q.value),]  # in case a gene presented in > 1 adjacent regions, select only the most sig. region
  }
  
  if (my_df$estimate > 0){
    return("increase")
  } else if (my_df$estimate < 0){
    return("decrease")
  }
}

direction <- lapply(1:nrow(gene_df), get_direction, gene_df = gene_df, gain = gain, loss = loss)
direction <- unlist(direction)

gene_df$direction <- direction
gene_df <- gene_df[,c("cancer_type", "gene", "GainOrLoss", "direction")]
gene_df$condition <- paste(gene_df$GainOrLoss, "-", gene_df$direction)
gene_df$GainOrLoss <- NULL
gene_df$direction <- NULL


df <- reshape2::dcast(gene_df, cancer_type ~ gene)
df[is.na(df)] <- ""


cols <- c("loss - decrease" = "#2b83ba", "loss - increase" = "#abdda4", 
          "gain - decrease" = "#fdae61", "gain - increase" = "#d7191c")

### Fig. 3d
pdf("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Age-associated_SCNA_in_cancer_gene.pdf", width = 12.5, height = 4) 
p <- ggplot(data = gene_df, aes(x=gene, y=cancer_type, fill=condition)) +
  geom_point(size = 3, colour="black", pch=21) +
  scale_fill_manual(values = cols, labels=c("Decrease gain", "Increase gain", 
                                            "Decrease loss", "Increase loss")) +
  scale_y_discrete(limits = rev(levels(gene_df$cancer_type))) +
  ggtitle("Age-associated SCNA changes in cancer driver genes") +
  xlab("Gene") + ylab("Cancer type") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9"))
print(p)
dev.off()

