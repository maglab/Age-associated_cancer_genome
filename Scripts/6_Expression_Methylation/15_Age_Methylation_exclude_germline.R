### Methylation with age by linear regression excluding samples with germline mutations

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)
library(biomaRt)

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")

# patients with germline variants
germline_patients <- read.table("Data/TCGA_PANCAN_Germline/TCGA_Patients_with_Germline_Variants.txt", header = FALSE)
germline_patients <- as.character(unlist(germline_patients))

# biomaRt
ensembl_100 <- useMart(host='http://apr2020.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
bmIDs <- getBM(attributes=c('external_gene_name', 'gene_biotype', 'entrezgene_id'), mart = ensembl_100)
bmIDs <- bmIDs[bmIDs$gene_biotype == "protein_coding",]

bmIDs <- bmIDs[complete.cases(bmIDs),]

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("BRCA")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("OV", "UCEC")){
    model <- gene ~ age + purity + race + figo_stage + histologic_grade
  } 
  return(model)
}


#### function to test linear regression for methylation of each gene
age_gene_methylation <- function(gene, project, df, clin, model){
  print(gene)
  df_tmp <- df[rownames(df) == gene,]
  df_tmp <- as.data.frame(t(df_tmp))
  colnames(df_tmp) <- "gene"
  df_tmp <- merge(df_tmp, clin, by.x = "row.names", by.y = "patient")
  colnames(df_tmp) <- c("patient", colnames(df_tmp)[2:ncol(df_tmp)])
  
  # model
  model <- as.formula(model)
  
  lm_fit <- lm(formula = model, data=df_tmp)  # fit linear model
  summary(lm_fit)
  result <- tidy(lm_fit)
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- gene
  colnames(result_df) <- c("gene", "estimate", "std.error", "statistic", "p.value")
  return(result_df)
}



### function to get age-DGEs from cancer
age_methylation_project <- function(project, germline_patients, purity){
  
  print(paste0("Working on: ", project))
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  clin <- merge(clin, purity, by.x = "patient", by.y = "patient")
  
  ### read expression result
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv"))
  genes_expression <- as.character(expression$gene)
  
  ### read methylation data
  df <- read.csv(paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", collapse = ""))
  print("Finish reading methylation table")
  
  rownames(df) <- df$X
  df$X <- NULL
  
  df <- df[rownames(df) %in% genes_expression,]
  genes <- rownames(df)
  
  colnames(df) <- unlist(lapply(colnames(df), clean_id))
  
  ### remove patients with germline variants
  df <- df[,!(colnames(df) %in% germline_patients)]
  dim(df)
  
  ### test linear regression for each gene
  print("Start linear regression")
  model <- format(terms(model_selection(project)))
  results <- lapply(genes, age_gene_methylation, project = project, df = df, clin = clin, model = model)
  results <- do.call(rbind, results)
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.csv(results, paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age_exclude_germline.csv", collapse = ""), row.names = FALSE)
  return(dim(df))
}

lapply(c("BRCA", "OV", "UCEC"), age_methylation_project, germline_patients = germline_patients, purity = purity)

