### get POLE and POLD1 Mutations across cancers from maf files
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(dplyr)

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


### function to generate table for samples and gene mutations
pol_mut_sample_table <- function(project){
  
  print(paste("Working on: ", project))
  
  # read maf file
  maf <- read.table(paste0("./Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""), header = TRUE)
  
  maf <- maf[maf$Short_ID %in% ASCAT_Tumour,] # select only samples in ASCAT
  num_samples <- length(unique(as.character(maf$Patient_ID)))
  
  all_samples <- unique(as.character(maf$Patient_ID))
  
  # filter to keep only variant we want
  maf <- maf[maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                               "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 
                                               "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"),]
  
  # PLOE and POLD1

  genes_of_interest <- c("POLE", "POLD1")
  
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
  print(paste0(project, "; POLE: ", sum(df_3$POLE), ", POLD1: ", sum(df_3$POLD1)))
  
  write.csv(df_3, paste0("Analysis_results/Mutations/6_POLE_POLD1/", project, "_pol_mutations.csv", collapse = ""))
  
  return(paste(project, dim(df_3)))
}

lapply(projects, pol_mut_sample_table)
