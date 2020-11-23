### From each cancer maf file
### Generate table for samples and gene mutations (rows = samples, columns = genes, 0: not mutated, 1: mutated)
### Filter only mutations in gene bodies, excluding mutations in intergenic regions, introns, 5'UTRs, 5'Flanks, 3'UTRs, and 3'Flanks
### Filter only genes in cancer genes
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
mut_sample_table <- function(project){
  
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
  
  # remove MSI-H tumours from COAD, READ, STAD, UCEC
  if(project %in% c("COAD", "READ", "STAD", "UCEC")){
    # read MSI status file
    MSI <- read.csv(paste0("Data/MSI_Status/", project, "_MSI_Status.csv", collapse = ""))
    MSI_patients <- as.character(MSI[MSI$mononucleotide_and_dinucleotide_marker_panel_analysis_status == "MSI-H",]$bcr_patient_barcode)
    maf <- maf[!(maf$Patient_ID %in% MSI_patients),]
  }
  
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
  
  df <- maf[,c("Hugo_Symbol","Patient_ID")]
  df$value <- rep(1, nrow(df))
  df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
  df$Patient_ID <-  as.character(df$Patient_ID)
  
  df_2 <- df %>% 
    dplyr::group_by(Hugo_Symbol, Patient_ID) %>% 
    dplyr::summarise(Sum=sum(value, na.rm = T)) %>%
    tidyr::spread(Hugo_Symbol, Sum, fill=0) %>%
    as.data.frame()
  
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
  
  print((apply(df_3, 2, sum)/num_samples) >= 0.05)
  
  write.csv(df_3, paste0("Analysis_results/Mutations/1_table_mutations_samples/", project, "_mutations_filter_hypermutated.csv", collapse = ""))
  
  return(paste(project, dim(df_3)))
}

lapply(projects, mut_sample_table)

