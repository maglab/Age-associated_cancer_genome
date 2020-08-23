### Check num samples that are in ASCAT and clinical

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

remove_last_char <- function(id){
  return(substr(id, start = 1, stop = nchar(id) - 1))
}

### function to read methylation file and check num samples
check_num_meth <- function(project, clinical, ASCAT){
  
  ### clinical
  clin <- clinical[clinical$cancer_type == project,]
  
  ### samples in ASCAT
  ASCAT_project <- ASCAT[ASCAT$cancer_type == project, ]
  ASCAT_Tumour <- as.character(ASCAT_project$barcodeTumour)
  ASCAT_Tumour <- unlist(lapply(ASCAT_Tumour, clean_id))
  ASCAT_Tumour <- unlist(lapply(ASCAT_Tumour, remove_last_char))
  
  ### read methylation data
  meth_df <- read.delim(paste0("Data/TCGA_DNA_Methylation/gdac_", project, "_methylation/", 
                               project, ".meth.by_min_expr_corr.data.txt", collapse = ""))
  meth_df$X <- NULL
  meth_df$X.1 <- NULL
  meth_df$X.2 <- NULL
  meth_df <- meth_df[2:nrow(meth_df),]
  rownames(meth_df) <- meth_df$Hybridization.REF
  meth_df$Hybridization.REF <- NULL
  indx <- sapply(meth_df, is.factor)
  meth_df[indx] <- lapply(meth_df[indx], function(x) as.numeric(as.character(x)))
  
  colnames(meth_df) <- unlist(lapply(colnames(meth_df), clean_id))
  
  ### select only samples in ASCAT
  meth_df <- meth_df[,colnames(meth_df) %in% ASCAT_Tumour]
  colnames(meth_df) <- unlist(lapply(colnames(meth_df), clean_id_1))
  
  ### select only samples in clinical
  meth_df <- meth_df[,colnames(meth_df) %in% clin$patient]
  write.csv(meth_df, paste0("Data/TCGA_DNA_Methylation_Filtered/", project, "_methylation.csv", 
                            collapse = ""))
  return(paste0(project, ": ", ncol(meth_df), collapse = ""))
}

lapply(projects, check_num_meth, clinical = clinical, ASCAT = ASCAT)

