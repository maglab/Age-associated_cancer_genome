### Calculate rank-based normalisation SCNA scores
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

# function to clean sample ID
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", x = id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

rank_norm_SCNA <- function(project){
  # read data
  df <- read.csv(paste0("Analysis_results/CNAs/2_SCNA_scores/raw_SCNA_scores/", project, "_raw_SCNA_score.csv", collapse = ""))
  rownames(df) <- df$X
  df$X <- NULL
  rank_norm_focal <- rank(df$focal_score, ties.method = "average")/nrow(df)
  rank_norm_chrom <- rank(df$chrom_score, ties.method = "average")/nrow(df)
  rank_norm_arm <- rank(df$arm_score, ties.method = "average")/nrow(df)
  
  df_norm <- cbind(rank_norm_focal, rank_norm_chrom, rank_norm_arm)
  rownames(df_norm) <- unlist(lapply(rownames(df), clean_id))
  sum_score <- apply(df_norm, 1, sum)
  df_norm <- cbind(df_norm,sum_score)
  write.csv(df_norm, paste0("Analysis_results/CNAs/2_SCNA_scores/normalised_SCNA_scores/", project, "_normalised_SCNA_score.csv", collapse = ""))
  return(project)
}

lapply(TCGA_projects, rank_norm_SCNA)

