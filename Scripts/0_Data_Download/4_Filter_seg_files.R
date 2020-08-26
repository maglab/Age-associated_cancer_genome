### filter seg files to keep only samples that were shared with ASCAT files (in a samples_of_interest file)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

# Read files
samples_of_interest <- read.csv("Data/samples_in_ASCAT_and_seg.csv")

### Filter CNA seg files to keep only samples with ASCAT profile
seg_files <- list.files("Data/TCGA_segfiles", full.names = TRUE, pattern = "TCGA")

# create function to filter seg files
filter_seg <- function(seg_file){
  tmp <- unlist(strsplit(seg_file, split = "\\/"))[3]
  tmp <- unlist(strsplit(tmp, split = "-"))[2]
  project <- unlist(strsplit(tmp, split = "_"))[1]
  
  seg <- read.table(seg_file, header = TRUE)
  seg_tumour <- seg[seg$Sample %in% as.character(samples_of_interest$barcodeTumour),]
  seg_normal <- seg[seg$Sample %in% as.character(samples_of_interest$barcodeNormal),]
  
  write.table(seg_tumour, paste0("Data/TCGA_segfiles/Seg_file_Tumour/", project, "_noCNV_seg_tumour.txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(seg_normal, paste0("Data/TCGA_segfiles/Seg_file_Normal/", project, "_noCNV_seg_normal.txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  nT <- length(unique(seg_tumour$Sample))
  nN <- length(unique(seg_normal$Sample))
  
  return(paste0(project, ",  num Tumour: ", nT, " ; num Normal: ", nN))
}


lapply(seg_files, filter_seg)



