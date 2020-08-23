### Methylation with age by linear regression

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

bmIDs <- bmIDs[complete.cases(bmIDs),]

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
age_methylation_project <- function(project, purity){
  
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
  
  ### test linear regression for each gene
  print("Start linear regression")
  model <- format(terms(model_selection(project)))
  results <- lapply(genes, age_gene_methylation, project = project, df = df, clin = clin, model = model)
  results <- do.call(rbind, results)
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.csv(results, paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""), row.names = FALSE)
}

lapply(projects, age_methylation_project, purity = purity)


