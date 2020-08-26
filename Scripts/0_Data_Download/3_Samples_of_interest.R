### Find samples that are present in both ASCAT output files and seg files
# then filter them to keep only samples with age information
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

seg_files <- list.files("Data/TCGA_segfiles", full.names = TRUE, pattern = "TCGA")
samples_seg <- c()
for(i in seg_files){
  tmp <- unlist(strsplit(i, split = "\\/"))[3]
  tmp <- unlist(strsplit(tmp, split = "-"))[2]
  project <- unlist(strsplit(tmp, split = "_"))[1]
  print(project)
  
  seg <- read.table(i, header = TRUE)
  print(length(unique(as.character(seg$Sample))))
  samples_seg <- c(samples_seg, as.character(unique(seg$Sample)))
}

purity <- read.table("Data/TCGA.purity.txt", header = TRUE) # read samples with tumour purity
dim(purity) # 9861

ASCAT <- read.table("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt", header = TRUE) # read samples with ASCAT profile
dim(ASCAT) # 9873

ASCAT <- ASCAT[ASCAT$patient %in% purity$patient,]
dim(ASCAT) # 9861

### seg file samples that are also in ASCAT files
samples_seg_in_ASCAT <- samples_seg[samples_seg %in% c(as.character(ASCAT$barcodeTumour), as.character(ASCAT$barcodeNormal))]
length(samples_seg_in_ASCAT)

### ASCAT samples that are also in seg files
ASCAT_in_seg <- ASCAT[ASCAT$barcodeTumour %in% samples_seg_in_ASCAT,]
dim(ASCAT_in_seg)   # 9751 samples

# check if they also have age information
clinical <- read.csv("Data/all_clin_indexed.csv")
clin_age <- clinical[,c("submitter_id", "age_at_index")]
clin_age <- clin_age[complete.cases(clin_age),]

ASCAT_in_seg <- ASCAT_in_seg[ASCAT_in_seg$patient %in% clin_age$submitter_id,]
dim(ASCAT_in_seg)   # 9678 samples persisted

write.csv(ASCAT_in_seg, "Data/samples_in_ASCAT_and_seg.csv", row.names = FALSE)
