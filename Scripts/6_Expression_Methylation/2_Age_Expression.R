### Expression with age by linear regression

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")
library(ggplot2)
library(broom)
library(biomaRt)

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# remove projects that have samples < 100
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", 
                                       "READ", "THYM", "UCS", "UVM"))]

# biomaRt
ensembl_100 <- useMart(host='http://apr2020.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
bmIDs <- getBM(attributes=c('external_gene_name', 'gene_biotype', 'entrezgene_id'), mart = ensembl_100)
bmIDs <- bmIDs[bmIDs$gene_biotype == "protein_coding",]

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- gene ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- gene ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- gene ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- gene ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- gene ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- gene ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- gene ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- gene ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- gene ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- gene ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- gene ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- gene ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- gene ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- gene ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- gene ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- gene ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- gene ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- gene ~ age + purity + gender + pathologic_stage
  }
  return(model)
}


#### function to test linear regression for expression of each gene
age_gene_expression <- function(gene, project, df_agg_filtered, clin, model){
  print(gene)
  df_tmp <- df_agg_filtered[rownames(df_agg_filtered) == gene,]
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
age_DEGs_project <- function(project, purity, bmIDs){
  
  print(paste0("Working on: ", project))
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  clin <- merge(clin, purity, by.x = "patient", by.y = "patient")
  
  ### read methylation data
  df_methylation <- read.csv(paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", collapse = ""))
  genes_methylation <- as.character(df_methylation$X)
  
  ### read expression data
  df <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered.csv", collapse = ""))
  print("Finish reading expression table")
  
  df <- df[df$gene %in% bmIDs$external_gene_name,]  # only protein-coding genes in biomaRt
  df <- df[df$gene %in% genes_methylation,]         # only genes in methylation
  
  # average duplicated genes
  df_agg <- aggregate(df[, -c(1)],
                      by = list(gene = df$gene),
                      FUN = mean,
                      na.rm = TRUE)
  print("Finish aggregate")
  dim(df_agg)
  rownames(df_agg) <- df_agg$gene
  df_agg$gene <- NULL
  
  keep.exprs <- rowSums(df_agg > 0) >= (0.5 * ncol(df_agg)) # filter out low count genes by keeping only genes with normalised RSEM > 1 in more than 30% of samples
  df_agg_filtered <- df_agg[keep.exprs,]
  dim(df_agg_filtered)
  df_agg_filtered <- log2(df_agg_filtered + 1)
  colnames(df_agg_filtered) <- unlist(lapply(colnames(df_agg_filtered), clean_id))
  genes <- rownames(df_agg_filtered)
  
  ### test linear regression for each gene
  print("Start linear regression")
  model <- format(terms(model_selection(project)))
  results <- lapply(genes, age_gene_expression, project = project, df_agg_filtered = df_agg_filtered, clin = clin, model = model)
  results <- do.call(rbind, results)
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.csv(results, paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv"), row.names = FALSE)
}

lapply(projects[!(projects %in% c("COAD", "READ", "UCEC"))], age_DEGs_project, purity = purity, bmIDs = bmIDs) # because COAD, READ, UCEC also have RNA-seq from Illumina GA


### Gene expression with age in COAD, READ, UCEC ##################################################

model_selection_1 <- function(project){
  if(project %in% c("COAD", "READ")){
    model <- gene ~ age + purity + gender + pathologic_stage + subtype + RNA_seq_platform
  } else if(project == "UCEC"){
    model <- gene ~ age + purity + race + figo_stage + histologic_grade + RNA_seq_platform
  }
  return(model)
}

age_gene_expression_1 <- function(gene, project, df_agg_filtered, clin, model, Hiseq_patients, GA_patients){
  print(gene)
  df_tmp <- df_agg_filtered[rownames(df_agg_filtered) == gene,]
  df_tmp <- as.data.frame(t(df_tmp))
  colnames(df_tmp) <- "gene"
  df_tmp <- merge(df_tmp, clin, by.x = "row.names", by.y = "patient")
  colnames(df_tmp) <- c("patient", colnames(df_tmp)[2:ncol(df_tmp)])
  df_tmp$RNA_seq_platform <- ifelse(df_tmp$patient %in% Hiseq_patients, "Hiseq", "GA")
  
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
age_DEGs_project <- function(project, purity, bmIDs){
  
  print(paste0("Working on: ", project))
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  clin <- merge(clin, purity, by.x = "patient", by.y = "patient")
  
  ### read methylation data
  df_methylation <- read.csv(paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", collapse = ""))
  genes_methylation <- as.character(df_methylation$X)
  
  ### read expression data
  df <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered.csv", collapse = ""))
  df_GA <- read.csv(paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered_GA_not_overlap_with_Hiseq.csv", collapse = ""))
  
  print("Finish reading expression table")
  
  df <- df[df$gene %in% bmIDs$external_gene_name,]  # only protein-coding genes in biomaRt
  df <- df[df$gene %in% genes_methylation,]         # only genes in methylation
  
  df_GA <- df_GA[df_GA$gene %in% bmIDs$external_gene_name,]  # only protein-coding genes in biomaRt
  df_GA <- df_GA[df_GA$gene %in% genes_methylation,]         # only genes in methylation
  
  # average duplicated genes
  df_agg <- aggregate(df[, -c(1)],
                      by = list(gene = df$gene),
                      FUN = mean,
                      na.rm = TRUE)
  df_GA_agg <- aggregate(df_GA[, -c(1)],
                      by = list(gene = df_GA$gene),
                      FUN = mean,
                      na.rm = TRUE)
  print("Finish aggregate")
  dim(df_agg)
  dim(df_GA_agg)
  rownames(df_agg) <- df_agg$gene
  rownames(df_GA_agg) <- df_GA_agg$gene
  df_agg$gene <- NULL
  df_GA_agg$gene <- NULL
  colnames(df_agg) <- unlist(lapply(colnames(df_agg), clean_id))
  colnames(df_GA_agg) <- unlist(lapply(colnames(df_GA_agg), clean_id))
  
  Hiseq_patients <- colnames(df_agg)
  GA_patients <- colnames(df_GA_agg)
  
  # merge both df
  df_agg <- cbind(df_agg, df_GA_agg)
  
  keep.exprs <- rowSums(df_agg > 0) >= (0.5 * ncol(df_agg)) # filter out low count genes by keeping only genes with normalised RSEM > 1 in more than 30% of samples
  df_agg_filtered <- df_agg[keep.exprs,]
  dim(df_agg_filtered)
  df_agg_filtered <- log2(df_agg_filtered + 1)
  colnames(df_agg_filtered) <- unlist(lapply(colnames(df_agg_filtered), clean_id))
  genes <- rownames(df_agg_filtered)
  
  ### test linear regression for each gene
  print("Start linear regression")
  model <- format(terms(model_selection_1(project)))
  
  results <- lapply(genes, age_gene_expression_1, project = project, df_agg_filtered = df_agg_filtered, 
                    clin = clin, model = model, Hiseq_patients = Hiseq_patients, GA_patients = GA_patients)
  results <- do.call(rbind, results)
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.csv(results, paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv"), row.names = FALSE)
}

lapply(c("COAD", "READ", "UCEC"), age_DEGs_project, purity = purity, bmIDs = bmIDs)





