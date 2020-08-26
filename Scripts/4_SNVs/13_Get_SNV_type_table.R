### Generate table for samples and SNV types (prepare table for oncoprint)
### Filter only mutations in gene bodies, excluding mutations in intergenic regions, introns, 5'UTRs, 5'Flanks, 3'UTRs, and 3'Flanks
### Filter only genes in cancer genes with 5% test
### Keep only non-silent mutations

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(dplyr)

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))
cancer_genes <- as.character(unlist(read.table("Data/cancer_genes.txt")))

# keep only projects that have > 100 samples
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "LAML", "MESO", "THYM", "UCS", "UVM"))]

# samples in ASCAT
ASCAT <- read.delim("Data/ASCAT_TCGA_filtered/summary.ascatTCGA.penalty70_filtered.txt")
ASCAT_Tumour <- ASCAT$barcodeTumour

clean_id <- function(ID){
  tmp <- unlist(strsplit(ID, split = "-"))[1:4]
  return(paste0(tmp, collapse = "-"))
}

clean_id_1 <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

ASCAT_Tumour <- unlist(lapply(as.character(ASCAT_Tumour), clean_id))

### helper function to check if >= 5% of samples harbouring mutation in a gene
five_percent <- function(gene, num_samples, maf){
  #print(gene)
  maf_tmp <- maf[maf$Hugo_Symbol == gene,]
  samples <- length(unique(as.character(maf_tmp$Short_ID)))
  percentage <- samples/num_samples
  result <- percentage >= 0.05
  return(result)
}

### function to generate table for samples and gene mutations
mut_type_sample_table <- function(project){
  
  print(paste("Working on: ", project))
  
  # read maf file
  maf <- read.table(paste0("./Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""), header = TRUE)
  
  maf <- maf[maf$Short_ID %in% ASCAT_Tumour,] # select only samples in ASCAT
  
  # read mut burden
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id_1))
  patient_ID_keep <- mut_burden_df[mut_burden_df$total < 1000, ]$Tumor_Sample_Barcode
  
  # filter only samples with less than 1000 mutations
  maf <- maf[maf$Patient_ID %in% patient_ID_keep,]
  
  num_samples <- length(unique(as.character(maf$Patient_ID)))
  
  all_samples <- unique(as.character(maf$Patient_ID))
  
  # filter to keep only variant we want
  maf <- maf[maf$Hugo_Symbol %in% cancer_genes,]
  maf <- maf[maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                               "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 
                                               "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"),]
  
  # genes which are mutated in >= 5% of samples
  genes <- unique(as.character(maf$Hugo_Symbol))
  five_percent_test <- lapply(genes, five_percent, num_samples = num_samples, maf = maf) # call function for 5% test
  names(five_percent_test) <- genes
  five_percent_result <- as.data.frame(do.call(rbind, five_percent_test))
  colnames(five_percent_result) <- "Pass"
  five_percent_result$Gene <- row.names(five_percent_result)
  row.names(five_percent_result) <- NULL
  five_percent_result <- five_percent_result[,c(2,1)]
  
  genes_of_interest <- five_percent_result[five_percent_result$Pass == TRUE,]$Gene
  
  print(paste("Genes of interest: ", paste0(genes_of_interest, collapse = " ")))
  
  # filter to keep only these genes
  maf <- maf[maf$Hugo_Symbol %in% genes_of_interest,]
  dim(maf)
  
  df <- maf[,c("Hugo_Symbol","Patient_ID", "Variant_Classification")]
  df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
  df$Patient_ID <- as.character(df$Patient_ID)
  df$Variant_Classification <- as.character(df$Variant_Classification)
  
  # remove duplicated SNVs (same type of SNVs occur in a gene in a patient)
  df_remove_dup <- df[!duplicated(df[c("Hugo_Symbol","Patient_ID", "Variant_Classification")]),]
  
  # check duplicated (the sample contains more than 1 type of SNVs in a gene) and convert them to "Multiple"
  dupe <- df_remove_dup[,c("Hugo_Symbol","Patient_ID")] # select columns to check duplicates
  df_remove_dup[rownames(df_remove_dup[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),]),"Variant_Classification"] <- "Multiple"
  df_remove_dup <- df_remove_dup[!duplicated(df_remove_dup[c("Hugo_Symbol","Patient_ID", "Variant_Classification")]),]
  
  df_2 <- reshape2::dcast(df_remove_dup, Hugo_Symbol ~ Patient_ID)
  df_2[is.na(df_2)] <- ""
  
  # add samples that do not have any mutation in our genes of interest
  samples_to_be_added <- all_samples[!(all_samples %in% colnames(df_2))]
  m <- matrix("", nrow = nrow(df_2), ncol = length(samples_to_be_added))
  m <- data.frame(m)
  colnames(m) <- samples_to_be_added
  rownames(m) <- rownames(df_2)
  
  df_3 <- cbind(df_2, m)
  
  write.csv(df_3, paste0("Analysis_results/Mutations/3.2_Age_SNVs_mut_types/", project, "_mutation_types_new.csv", collapse = ""), row.names = FALSE)
  
  return(paste(project, dim(df_3)))
}

lapply(projects, mut_type_sample_table)


