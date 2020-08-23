### Correlation between copy-number from GISTIC2 and gene expression
### regression by dividing samples into 1) "Highly Loss", 2) "Loss", 3) "No Change", 4) "Gain", 5) "Highly Gain"
### then modelled gene expression using CNA status and purity
### Fig. 3e
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(data.table)
library(broom)


# read genes of interest
gain <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_gain.csv")
loss <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_del.csv")

cancer_type_gain <- sort(unique(as.character(gain$cancer_type)))
cancer_type_loss <- sort(unique(as.character(loss$cancer_type)))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:4]
  return(paste0(tmp, collapse = "-"))
}

clean_id1 <- function(id){
  tmp <- unlist(strsplit(id, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}



### function to test correlation for each gene
gene_corr <- function(gene, cancer_type, GISTIC, expression, GainOrLoss){
  print(paste0("working on: ", gene))
  tmp_expression <- as.data.frame(t(expression[rownames(expression) == gene,]))
  row.names(tmp_expression) <- lapply(row.names(tmp_expression), gsub, pattern = "[.]", replacement = "-")
  colnames(tmp_expression) <- "expression"
  
  tmp_GISTIC <- as.data.frame(t(GISTIC[rownames(GISTIC) == gene,]))
  colnames(tmp_GISTIC) <- "GISTIC"
  
  tmp_df <- merge(tmp_GISTIC, tmp_expression, by = "row.names")
  colnames(tmp_df) <- c("patient", "GISTIC", "expression")
  
  corr <- cor.test(tmp_df$GISTIC, tmp_df$expression, method = "pearson")
  
  p_value <- as.numeric(corr$p.value)
  p_plot <- formatC(as.numeric(corr$p.value), format = "e", digits = 2)
  
  r <- corr$estimate
  r_plot <- round(as.numeric(corr$estimate), 2)
  
  ### GISTIC cutoff
  tmp_df$GISTIC[tmp_df$GISTIC >= 1] <- 2
  tmp_df$GISTIC[tmp_df$GISTIC < 1 & tmp_df$GISTIC >= 0.25] <- 1
  tmp_df$GISTIC[tmp_df$GISTIC >= -0.25 & tmp_df$GISTIC < 1] <- 0
  tmp_df$GISTIC[tmp_df$GISTIC >= -1 & tmp_df$GISTIC < -0.25] <- -1
  tmp_df$GISTIC[tmp_df$GISTIC < -1] <- -2
  
  ### prepare table for plotting
  new_labels <- c()
  for(i in 1:length(tmp_df$GISTIC)){
    if(tmp_df$GISTIC[i] == -2){
      new_labels <- c(new_labels, "Highly Loss")
    } else if(tmp_df$GISTIC[i] == -1){
      new_labels <- c(new_labels, "Loss")
    } else if(tmp_df$GISTIC[i] == 0){
      new_labels <- c(new_labels, "No Change")
    } else if(tmp_df$GISTIC[i] == 1){
      new_labels <- c(new_labels, "Gain")
    } else if(tmp_df$GISTIC[i] == 2){
      new_labels <- c(new_labels, "Highly Gain")
    }
  }
  tmp_df$CNA_status <- new_labels
  tmp_df$CNA_status <- as.factor(tmp_df$CNA_status)
  all_levels <- c("Highly Loss", "Loss", "No Change", "Gain", "Highly Gain")
  missing_levels <- all_levels[!(all_levels %in% tmp_df$CNA_status)]
  levels(tmp_df$CNA_status) <- c(levels(tmp_df$CNA_status), missing_levels)
  tmp_df$CNA_status <- ordered(tmp_df$CNA_status, levels = all_levels)  # reorder levels
  
  # plot
  my_label <- paste0("r = ", r_plot, "\np = ", p_plot)
  p <- ggplot(tmp_df, aes(x = CNA_status, y = expression, fill = CNA_status, group = CNA_status)) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(width = 0.2, fill = "white") +
    scale_fill_brewer(palette="RdYlBu", drop = FALSE, direction = -1) + 
    scale_x_discrete(limits = c("Highly Loss", "Loss", "No Change", "Gain", "Highly Gain")) +
    labs(title = paste0(cancer_type, ": ", gene), 
         x = "SCNAs", y = "log2(normalised RSEM + 1)") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size=11,face="bold"),
          axis.title.y = element_text(size=11,face="bold"),
          legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 5,
             label = my_label, hjust=0, vjust=1)
  
  pdf(paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/SCNAs_gene_expression/", 
             cancer_type, "_", gene, "_", GainOrLoss, "_expression.pdf", collapse = ""), 
      width = 6, height = 4.5) 
  print(p)
  dev.off()
  
  return(tidy(corr))
}

### function to test every gene for each cancer
gene_corr_cancer <- function(cancer_type, GainOrLoss, df){
  
  print(paste0("Working on cancer type: ", cancer_type))
  ### select only cancer of interest
  my_df <- df[df$cancer_type == cancer_type,]
  my_df$genes
  
  genes <- c()
  for(i in 1:length(my_df$genes)){
    region_genes <- my_df$genes[i]
    genes <- c(genes, unlist(strsplit(as.character(region_genes), split = ",")))
  }
  genes <- unique(genes)
  
  ### read GISTIC result files to get sample names
  GISTIC <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", cancer_type,
                              "_GISTIC_mkCNV/", "all_data_by_genes.txt", collapse = ""), sep = "\t")
  
  GISTIC <- GISTIC[GISTIC$Gene.Symbol %in% genes,]
  row.names(GISTIC) <- GISTIC$Gene.Symbol
  GISTIC$Gene.Symbol <- NULL
  GISTIC$Gene.ID <- NULL
  GISTIC$Cytoband <- NULL
  
  colnames(GISTIC) <- unlist(lapply(colnames(GISTIC), clean_id))
  
  ### read expression table
  expression <- read.csv(paste0("Data/TCGA_Gene_Expression/TCGA-", cancer_type, 
                                "_normalised_expression.csv", collapse = ""))
  
  ### clean ID in expression table
  IDs <- colnames(expression)[2:ncol(expression)]
  IDs <- unlist(lapply(IDs, clean_id))
  colnames(expression) <- c("gene", IDs)
  
  expression <- expression[expression$gene %in% genes, ]
  expression <- data.table(expression)
  
  ### average expression if a gene is presented in > 1 rows
  # average duplicated rows (if a gene is presented in > 1 rows)
  expression <- expression[ ,lapply(.SD, mean), by = gene]
  expression <- as.data.frame(expression)
  
  rownames(expression) <- expression$gene
  expression$gene <- NULL
  
  ### samples shared by both GISTIC and gene expression
  shared_IDs <- IDs[IDs %in% colnames(GISTIC)]
  GISTIC <- GISTIC[,colnames(GISTIC) %in% shared_IDs]
  expression <- expression[, colnames(expression) %in% shared_IDs]
  
  ### genes shared by both GISTIC and gene expression
  shared_genes <- rownames(GISTIC)[rownames(GISTIC) %in% rownames(expression)]
  GISTIC <- GISTIC[rownames(GISTIC) %in% shared_genes,]
  expression <- expression[rownames(expression) %in% shared_genes,]
  
  ### clean id
  colnames(GISTIC) <- unlist(lapply(colnames(GISTIC), clean_id1))
  colnames(expression) <- unlist(lapply(colnames(expression), clean_id1))
  
  ### log transform expression
  expression <- log2(expression + 1)
  
  ### call the function to test every genes
  results <- lapply(shared_genes, gene_corr, cancer_type = cancer_type, expression = expression, 
                    GISTIC = GISTIC, GainOrLoss = GainOrLoss)
  names(results) <- shared_genes
  
  results <- as.data.frame(do.call(rbind,results))
  results$adj_p <- p.adjust(results$p.value, method = "BH")
  results$sig <- ifelse(results$adj_p < 0.05 & results$estimate > 0, TRUE, FALSE)
  
  write.csv(results, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/SCNAs_gene_expression/", 
                            cancer_type, "_CNAs_Gene_Expression_", GainOrLoss, ".csv", collapse = ""), row.names = TRUE)
}


# for gain
lapply(cancer_type_gain, gene_corr_cancer, GainOrLoss = "gain", df = gain)

# for loss
lapply(cancer_type_loss, gene_corr_cancer, GainOrLoss = "loss", df = loss)

