### Keep only samples in ASCAT output and those that have age information

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")


### Read data
# clinical
clinical <- read.csv("Data/all_clin_XML.csv")

# samples in ASCAT
ASCAT <- read.delim("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt")
ASCAT_Tumour <- as.character(ASCAT$barcodeTumour)

# clean id 
clean_id <- function(ID){
  tmp <- unlist(strsplit(ID, split = "-"))[1:4]
  tmp[4] <- substr(tmp[4], 1, 2)
  return(paste0(tmp, collapse = "-"))
}

clean_id_1 <- function(id){
  tmp <- unlist(strsplit(id, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

ASCAT_Tumour <- unlist(lapply(ASCAT_Tumour, clean_id))

### Pathway alteration data
pathway_df <- read.csv("Data/Pathway_alterations.csv")

pathway_df <- pathway_df[as.character(pathway_df$SAMPLE_BARCODE) %in% ASCAT_Tumour,]
dim(pathway_df)

pathway_df$SAMPLE_BARCODE <- unlist(lapply(as.character(pathway_df$SAMPLE_BARCODE), clean_id_1))
pathway_df <- pathway_df[pathway_df$SAMPLE_BARCODE %in% clinical$patient,]
dim(pathway_df)  # 8055 samples

### Check num samples per cancer type
clinical_tmp <- clinical[clinical$patient %in% pathway_df$SAMPLE_BARCODE,]
table(clinical_tmp$cancer_type)

projects <- names(table(clinical_tmp$cancer_type)[table(clinical_tmp$cancer_type) > 100])

write.csv(pathway_df, "Data/Pathway_alterations_clean.csv", row.names = FALSE)
