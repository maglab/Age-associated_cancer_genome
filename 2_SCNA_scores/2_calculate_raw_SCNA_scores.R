### Calculate SCNA score (Method from Yuan et al., Cancer Cell, 2018)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

### Focal-Level ###################################################################################
# function to return cutoff score (-2, -1, 0, 1, 2) from raw copy numbers
get_score <- function(score){
  if(score < -1){
    return(-2)
  } else if(score < -0.25 & score >= -1){
    return(-1)
  } else if(score >= -0.25 & score < 0.25){
    return(0)
  } else if(score >= 0.25 & score < 1){
    return(1)
  } else if(score >= 1){
    return(2)
  }
}

# function to return the score for each patient (sum of the absolute values of every focal region)
sum_score <- function(score_vector){
  return(sum(unlist(lapply(score_vector, abs))))
}

# wrap up get_score and sum_score
score_for_patient <- function(raw_score_vector){
  return(sum_score(unlist(lapply(raw_score_vector, get_score))))
}

### Main function to return focal-score for all patients in a project
focal_score <- function(project){
  # read GISTIC focal-level result file
  df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt"))
  df$X <- NULL
  df <- df[df$Amplitude.Threshold == "Actual Copy Change Given",]   # keep only raw values
  df <- df[,-c(1:9)]    # remove unnecessary columns
  df[1:10,1:10]
  
  result <- apply(df, 2, score_for_patient)
  
  return(result)
}

### Arm- and Chromosome-Level #####################################################################
# function to check which arms should be considered as a chromo-level score
check_chrom <- function(score, arms){
  tmp <- as.data.frame(cbind(arms, score = score))
  tmp$chrom <- unlist(lapply(as.character(tmp$arms), gsub, pattern = "[pq]", replacement = ""))
  
  chrom_level <- c()  # vector to store which arms should be considered as a chrom-level
  for(i in unique(tmp$chrom)){
    
    tmp_chrom <- tmp[tmp$chrom == i,]
    
    if(length(tmp_chrom$score) == 2){
      if(tmp_chrom$score[1] == tmp_chrom$score[2] & tmp_chrom$score[1] != 0){
        chrom_level <- c(chrom_level, as.character(tmp_chrom$arms[1]))
        chrom_level <- c(chrom_level, as.character(tmp_chrom$arms[2]))
      }
    }
  }
  
  if(length(chrom_level) == 0){
    return("No_chrom")
  } else{
    return(chrom_level)
  }
}

# function to calculate score for chrom-level
chrom_score <- function(for_chrom, score_df){
  tmp <- score_df[score_df$arms %in% for_chrom, ]
  tmp <- tmp[duplicated(tmp$chrom),]
  score_for_chrom <- unlist(lapply(as.numeric(as.character(tmp$score)), get_score))
  return(sum_score(score_for_chrom))
}

# function to calculate score for arm-level
arm_score <- function(for_arm, score_df){
  tmp <- score_df[score_df$arms %in% for_arm, ]
  score_for_arm <- unlist(lapply(as.numeric(as.character(tmp$score)), get_score))
  return(sum_score(score_for_arm))
}

# wrap up function for each patient
patient_chrom_arm <- function(score, arms){
  for_chrom <- check_chrom(score, arms)   # arms considered in chrom-level
  for_arm <- arms[!(arms %in% for_chrom)] # arms considered in arm-level
  
  score_df <- as.data.frame(cbind(arms, score = score))
  score_df$chrom <- unlist(lapply(as.character(score_df$arms), gsub, pattern = "[pq]", replacement = ""))
  
  if("No_chrom" %in% for_chrom){
    score_for_chrom <- 0
    score_for_arm <- arm_score(for_arm, score_df)
  } else {
    score_for_chrom <- chrom_score(for_chrom, score_df)
    score_for_arm <- arm_score(for_arm, score_df)
  }
  return(cbind(chrom_score = score_for_chrom, arm_score = score_for_arm))
}

### Main function to return chrom- and arm-score for all patients in a project
chrom_arm_score <- function(project){
  # read GISTIC arm-level result file
  df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/broad_values_by_arm.txt"))
  rownames(df) <- df$Chromosome.Arm
  df$Chromosome.Arm <- NULL
  df[1:10,1:10]
  
  arms <- rownames(df)
  
  result <- apply(df, 2, patient_chrom_arm, arms = arms)
  result <- t(result)
  colnames(result) <- c("chrom_score", "arm_score")
  return(result)
}

### Join Focal, Chrom, Arm scores together ########################################################
SCNA_raw_score <- function(project){
  focal <- focal_score(project)
  chrom_arm <- chrom_arm_score(project)
  my_df <- cbind(focal_score = focal, chrom_arm)
  write.csv(my_df, paste0("Analysis_results/CNAs/2_SCNA_scores/raw_SCNA_scores/", project, "_raw_SCNA_score.csv", collapse = ""))
}

lapply(TCGA_projects, SCNA_raw_score)

