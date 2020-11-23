### Mutation table for PANCAN analysis

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(dplyr)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))
cancer_genes <- as.character(unlist(read.table("Data/cancer_genes.txt")))

# read mc3 maf file
maf_file <- read.table("Data/mc3_with_patient_barcode.maf", header = TRUE)

# samples in ASCAT
ASCAT <- read.delim("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt")
ASCAT_Tumour <- ASCAT$barcodeTumour

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:4]
  return(paste0(tmp, collapse = "-"))
}

clean_id_1 <- function(ID){
  tmp <- unlist(strsplit(ID, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

ASCAT_Tumour <- unlist(lapply(as.character(ASCAT_Tumour), clean_id))

### select only samples in ASCAT
maf_file$Short_ID <- unlist(lapply(as.character(maf_file$Tumor_Sample_Barcode), clean_id))
maf_file_1 <- maf_file[maf_file$Short_ID %in% ASCAT_Tumour,]

### select only samples in clinical file
ids <- unique(as.character(clinical$patient))

maf_file_1 <- maf_file_1[maf_file_1$Patient_ID %in% ids,]

### function to read mut burden for each cancer type and combine together
get_mut_burden <- function(project){
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id_1))
  mut_burden_df <- mut_burden_df[,c("Tumor_Sample_Barcode", "total")]
  mut_burden_df <- mut_burden_df[mut_burden_df$total < 1000,]
  print(paste0(project, ": ", nrow(mut_burden_df)))
  return(mut_burden_df)
}

PANCAN_mut <- lapply(projects, get_mut_burden)
PANCAN_mut <- do.call(rbind, PANCAN_mut)
dim(PANCAN_mut) # 8585 samples

# remove MSI-H tumours from COAD, READ, STAD, UCEC
MSI_COAD <- read.csv("Data/MSI_Status/COAD_MSI_Status.csv")
MSI_COAD <- as.character(MSI_COAD[MSI_COAD$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)
MSI_READ <- read.csv("Data/MSI_Status/READ_MSI_Status.csv")
MSI_READ <- as.character(MSI_READ[MSI_READ$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)
MSI_STAD <- read.csv("Data/MSI_Status/STAD_MSI_Status.csv")
MSI_STAD <- as.character(MSI_STAD[MSI_STAD$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)
MSI_UCEC <- read.csv("Data/MSI_Status/UCEC_MSI_Status.csv")
MSI_UCEC <- as.character(MSI_UCEC[MSI_UCEC$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)

MSI_patients <- c(MSI_COAD, MSI_READ, MSI_STAD, MSI_UCEC)

PANCAN_mut <- PANCAN_mut[!(PANCAN_mut$Tumor_Sample_Barcode %in% MSI_patients),]
dim(PANCAN_mut) # 8448 samples

### select only samples with mutation loads < 1000 and not MSI-H
maf_file_1 <- maf_file_1[maf_file_1$Patient_ID %in% PANCAN_mut$Tumor_Sample_Barcode,]
num_samples <- length(unique(as.character(maf_file_1$Patient_ID)))
num_samples

all_samples <- unique(as.character(maf_file_1$Patient_ID))

# filter to keep only variant we want
maf_file_1 <- maf_file_1[maf_file_1$Hugo_Symbol %in% cancer_genes,]
maf_file_1 <- maf_file_1[maf_file_1$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                                                  "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 
                                                                  "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"),]
dim(maf_file_1)


genes <- sort(unique(as.character(maf_file_1$Hugo_Symbol)))

### function to check if >= 5% of samples harbouring mutation in a gene
five_percent <- function(gene, num_samples, maf){
  print(gene)
  maf_tmp <- maf[maf$Hugo_Symbol == gene,]
  samples <- length(unique(as.character(maf_tmp$Short_ID)))
  percentage <- samples/num_samples
  result <- percentage >= 0.05
  return(result)
}

five_percent_result <- lapply(genes, five_percent, num_samples = num_samples, maf = maf_file_1)
names(five_percent_result) <- genes
genes_of_interest <- as.data.frame(do.call(rbind, five_percent_result))
genes_of_interest$gene <- row.names(genes_of_interest)
row.names(genes_of_interest) <- NULL
colnames(genes_of_interest) <- c("Pass", "gene")
genes_of_interest <- genes_of_interest[,c(2,1)]
head(genes_of_interest)
genes_of_interest <- genes_of_interest[genes_of_interest$Pass == TRUE,]$gene
length(genes_of_interest)  # 20 genes


# filter to keep only these genes
maf_file_2 <- maf_file_1[maf_file_1$Hugo_Symbol %in% genes_of_interest,]
dim(maf_file_2)

df <- maf_file_2[,c("Hugo_Symbol","Patient_ID")]
df$value <- rep(1, nrow(df))
df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
df$Patient_ID <-  as.character(df$Patient_ID)

df_2 <- df %>% 
  dplyr::group_by(Hugo_Symbol, Patient_ID) %>% 
  dplyr::summarise(Sum=sum(value, na.rm = T)) %>%
  tidyr::spread(Hugo_Symbol, Sum, fill=0) %>%
  as.data.frame()
head(df_2)

row.names(df_2) <- df_2$Patient_ID
df_2$Patient_ID <- NULL
df_2[df_2 > 0] <- 1   # in case some patients have mutations in more than 1 regions in a same gene

# add samples that do not have any mutation in our genes of interest
samples_to_be_added <- all_samples[!(all_samples %in% row.names(df_2))]
m <- matrix(0, ncol = ncol(df_2), nrow = length(samples_to_be_added))
m <- data.frame(m)
row.names(m) <- samples_to_be_added
colnames(m) <- colnames(df_2)

df_3 <- rbind(df_2, m)
dim(df_3)

write.csv(df_3, "Analysis_results/Mutations/1_table_mutations_samples/PANCAN_mutations.csv")
