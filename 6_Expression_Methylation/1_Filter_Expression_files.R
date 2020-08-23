### Check num samples that are in ASCAT output files & clinical

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# samples in ASCAT
ASCAT <- read.delim("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt")

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:4]
  return(paste0(tmp, collapse = "-"))
}

clean_id_1 <- function(id){
  tmp <- unlist(strsplit(id, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to read expression file and check num samples
check_num_expression <- function(project, clinical, ASCAT){
  
  print(paste0("Working on: ", project))
  ### samples in ASCAT
  ASCAT_project <- ASCAT[ASCAT$cancer_type == project, ]
  ASCAT_Tumour <- as.character(ASCAT_project$barcodeTumour)
  ASCAT_Tumour <- unlist(lapply(ASCAT_Tumour, clean_id))
  
  ### clinical
  clin <- clinical[clinical$cancer_type == project, ]
  
  ### read expression file
  df <- read.csv(paste0("Data/TCGA_Gene_Expression/TCGA-", project, "_normalised_expression.csv", collapse = ""))
  tmp_col <- colnames(df)
  tmp_col <- unlist(lapply(tmp_col, clean_id))
  colnames(df) <- c("gene", tmp_col[2:length(tmp_col)])
  
  keep <- tmp_col[tmp_col %in% ASCAT_Tumour]
  keep <- keep[unlist(lapply(keep, clean_id_1)) %in% clin$patient]
  df <- df[colnames(df) %in% c("gene", keep)]
  
  ### remove duplicated samples to keep only one (eg. remove "TCGA-21-1076-01A.1" and keep "TCGA-21-1076-01A" from LUSC)
  df <- df[,!grepl("[.]", colnames(df))]
  
  write.csv(df, paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered.csv", collapse = ""), row.names = FALSE)
  return(paste0(project, ": ", ncol(df)-1))
}

lapply(projects, check_num_expression, clinical = clinical, ASCAT = ASCAT)


#################################### For Illumina GA platform #####################################
### function to read expression file and check num samples
check_num_expression_GA <- function(project, clinical, ASCAT){
  
  print(paste0("Working on: ", project))
  ### samples in ASCAT
  ASCAT_project <- ASCAT[ASCAT$cancer_type == project, ]
  ASCAT_Tumour <- as.character(ASCAT_project$barcodeTumour)
  ASCAT_Tumour <- unlist(lapply(ASCAT_Tumour, clean_id))
  
  ### clinical
  clin <- clinical[clinical$cancer_type == project, ]
  
  ### read expression file
  df <- read.csv(paste0("Data/TCGA_Gene_Expression/TCGA-", project, "_normalised_expression_GA.csv", collapse = ""))
  tmp_col <- colnames(df)
  tmp_col <- unlist(lapply(tmp_col, clean_id))
  colnames(df) <- c("gene", tmp_col[2:length(tmp_col)])
  
  keep <- tmp_col[tmp_col %in% ASCAT_Tumour]
  keep <- keep[unlist(lapply(keep, clean_id_1)) %in% clin$patient]
  df <- df[colnames(df) %in% c("gene", keep)]
  
  ### remove duplicated samples to keep only one (eg. remove "TCGA-21-1076-01A.1" and keep "TCGA-21-1076-01A" from LUSC)
  df <- df[,!grepl("[.]", colnames(df))]
  
  write.csv(df, paste0("Data/TCGA_Gene_Expression_Filtered/TCGA-", project, "_normalised_expression_filtered_GA.csv", collapse = ""), row.names = FALSE)
  return(paste0(project, ": ", ncol(df)-1))
}

lapply(c("COAD", "READ", "UCEC"), check_num_expression_GA, clinical = clinical, ASCAT = ASCAT)

### combine Illumina Hiseq and GA for COAD READ UCEC
clean_id_2 <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

COAD_Hiseq <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-COAD_normalised_expression_filtered.csv")
COAD_GA <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-COAD_normalised_expression_filtered_GA.csv")

COAD_GA <- COAD_GA[, c("gene", colnames(COAD_GA)[!(colnames(COAD_GA) %in% colnames(COAD_Hiseq))])]
colnames(COAD_GA) <- c("gene", unlist(lapply(colnames(COAD_GA[2:ncol(COAD_GA)]),clean_id_2)))
write.csv(COAD_GA, "Data/TCGA_Gene_Expression_Filtered/TCGA-COAD_normalised_expression_filtered_GA_not_overlap_with_Hiseq.csv", row.names = FALSE)


READ_Hiseq <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-READ_normalised_expression_filtered.csv")
READ_GA <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-READ_normalised_expression_filtered_GA.csv")

READ_GA <- READ_GA[, c("gene", colnames(READ_GA)[!(colnames(READ_GA) %in% colnames(READ_Hiseq))])]
colnames(READ_GA) <- c("gene", unlist(lapply(colnames(READ_GA[2:ncol(READ_GA)]),clean_id_2)))
write.csv(READ_GA, "Data/TCGA_Gene_Expression_Filtered/TCGA-READ_normalised_expression_filtered_GA_not_overlap_with_Hiseq.csv", row.names = FALSE)


UCEC_Hiseq <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-UCEC_normalised_expression_filtered.csv")
UCEC_GA <- read.csv("Data/TCGA_Gene_Expression_Filtered/TCGA-UCEC_normalised_expression_filtered_GA.csv")

UCEC_GA <- UCEC_GA[, c("gene", colnames(UCEC_GA)[!(colnames(UCEC_GA) %in% colnames(UCEC_Hiseq))])]
colnames(UCEC_GA) <- c("gene", unlist(lapply(colnames(UCEC_GA[2:ncol(UCEC_GA)]),clean_id_2)))
write.csv(UCEC_GA, "Data/TCGA_Gene_Expression_Filtered/TCGA-UCEC_normalised_expression_filtered_GA_not_overlap_with_Hiseq.csv", row.names = FALSE)


