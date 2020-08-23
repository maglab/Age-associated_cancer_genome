### script for separating mc3 maf files to multiple files based on TCGA projects

# maf file (mc3.v0.2.8.PUBLIC.maf) download from https://gdc.cancer.gov/about-data/publications/mc3-2017

# this will read mc3 maf file, add patient ID, and seperate it into several maf files by TCGA project
# keep only samples in clinical file (because they have age information) and in ASCAT output
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

# add patient ID to maf file
maf_file <- read.table("Data/mc3.v0.2.8.PUBLIC.maf", header = TRUE)

ids_tumour <- as.character(maf_file$Tumor_Sample_Barcode)

sample_id <- function(id){
  tmp <- unlist(strsplit(id, split = "-"))
  return(paste0(tmp[1:3], collapse = "-"))
}

maf_file$Patient_ID <- unlist(lapply(ids_tumour, sample_id))

write.table(maf_file, "Data/mc3_with_patient_barcode.maf", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

### separate by cancer type
# read mc3 maf file
maf_file <- read.table("Data/mc3_with_patient_barcode.maf", header = TRUE)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# samples in ASCAT
ASCAT <- read.delim("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt")
ASCAT_Tumour <- ASCAT$barcodeTumour

clean_id <- function(ID){
  tmp <- unlist(strsplit(ID, split = "-"))[1:4]
  return(paste0(tmp, collapse = "-"))
}

ASCAT_Tumour <- unlist(lapply(as.character(ASCAT_Tumour), clean_id))

### select only samples in ASCAT
maf_file$Short_ID <- unlist(lapply(as.character(maf_file$Tumor_Sample_Barcode), clean_id))
maf_file_1 <- maf_file[maf_file$Short_ID %in% ASCAT_Tumour,]

for(project in projects){
  ids <- unique(as.character(clinical[clinical$cancer_type == project,]$patient))
  maf_file_new <- maf_file_1[maf_file_1$Patient_ID %in% ids,]
  write.table(maf_file_new, paste0("./Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  print(c(project, length(unique(maf_file_new$Tumor_Sample_Barcode))))
}

# check num samples
get_num_samples <- function(project){
  maf <- read.table(paste0("./Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""), header = TRUE)
  return(length(as.character(unique(maf$Short_ID))))
}

lapply(projects, get_num_samples)

# check mutational burden
library(maftools)
get_mut_burden <- function(project){
  maf <- read.maf(paste0("Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""))
  maf_sample_summary <- getSampleSummary(maf)
  maf_sample_summary <- as.data.frame(maf_sample_summary)
  write.csv(maf_sample_summary, paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""), row.names = FALSE)
  return(paste(project, nrow(maf_sample_summary)))
}

lapply(projects,get_mut_burden)
