### get genes in the focal-regions significantly associated with age
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(biomaRt)

# read significant results
gain <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain_sig.csv")
del <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_del_sig.csv")

# get wide peak limits
wide_peak <- c()
for(i in 1:nrow(gain)){
  cancer_type <- gain$cancer_type[i]
  peak <- as.character(gain$Unique.Name[i])
  conf_region <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", cancer_type, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
  peak <- as.character(conf_region[as.character(conf_region$Unique.Name) == peak,]$Wide.Peak.Limits)
  peak <- gsub(pattern = " ", replacement = "", peak)
  wide_peak <- c(wide_peak, peak)
}

gain$WidePeak <- wide_peak


wide_peak <- c()
for(i in 1:nrow(del)){
  cancer_type <- del$cancer_type[i]
  peak <- as.character(del$Unique.Name[i])
  conf_region <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", cancer_type, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
  peak <- as.character(conf_region[as.character(conf_region$Unique.Name) == peak,]$Wide.Peak.Limits)
  peak <- gsub(pattern = " ", replacement = "", peak)
  wide_peak <- c(wide_peak, peak)
}

del$WidePeak <- wide_peak


gain_df <- gain[c("cancer_type", "region", "WidePeak")]
del_df <- del[c("cancer_type", "region", "WidePeak")]

# clean ID
clean_id <- function(region_name){
  return(unlist(strsplit(region_name, split = "X"))[2])
}

# get all protein-coding genes from biomaRt
# Select only protein coding genes using biomart (version 100, April 2020)
# biomart version 100 (April 2020)
ensembl <- useMart(host='http://apr2020.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL', 
                   dataset='hsapiens_gene_ensembl')
bmIDs <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype','entrezgene_id'),mart = ensembl)

protein_coding <- bmIDs[bmIDs$gene_biotype == 'protein_coding',]
protein_coding <- unique(protein_coding$external_gene_name)
length(protein_coding)  # all protein-coding genes

### function to extract gene from each region #####################################################
get_genes <- function(cancer_type, GainOrDel, region){
  # read GISTIC result of amp or del genes
  if(GainOrDel == "gain"){
    file <- "amp_genes.conf_95.txt"
  } else if(GainOrDel == "del"){
    file <- "del_genes.conf_95.txt"
  }
  gene_file <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", cancer_type, "_GISTIC_mkCNV/", file, collapse = ""))
  gene_file$X <- NULL
  colnames(gene_file) <- c(colnames(gene_file)[1], lapply(colnames(gene_file)[2:ncol(gene_file)], clean_id))
  gene_file <- gene_file[-c(1:3),] # remove unnecessary rows
  
  # get genes from region of interest
  genes <- as.character(unlist(gene_file[[region]]))  
  genes <- genes[genes != ""]
  genes <- unlist(lapply(genes, gsub, pattern = "\\[|\\]", replacement = ""))
  return(genes)
}

# function to keep only protein-coding genes from a gene list
keep_protein_coding <- function(genes){
  return(genes[genes %in% protein_coding])
}

# gain
gene_list <- mapply(get_genes, as.character(gain_df$cancer_type), GainOrDel = "gain", region = as.character(gain_df$region))
gene_list <- lapply(gene_list, keep_protein_coding)
gene_list <- lapply(gene_list, paste0, collapse = ",")
names(gene_list) <- NULL

gain <- cbind(gain, genes = do.call(rbind, gene_list))
View(gain)

# del
gene_list <- mapply(get_genes, as.character(del_df$cancer_type), GainOrDel = "del", region = as.character(del_df$region))
gene_list <- lapply(gene_list, keep_protein_coding)
gene_list <- lapply(gene_list, paste0, collapse = ",")
names(gene_list) <- NULL

del <- cbind(del, genes = do.call(rbind, gene_list))
View(del)

write.csv(gain, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_gain.csv", row.names = FALSE)
write.csv(del, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Genes_in_sig_focal_del.csv", row.names = FALSE)

