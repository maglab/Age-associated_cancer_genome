### Age biases in focal-level recurrent SCNAs (simple logistic regression)
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(broom)
library(GenomicRanges)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

gap <- read.table("Data/gap.txt")   # genome gap
gap <- gap[,c(2,3,4,7,8)]
colnames(gap) <- c("chr", "start", "end", "size", "type")
centro_telomere <- gap[gap$type %in% c("centromere", "telomere"),]   # select regions for centromere and telomere
dim(centro_telomere)

# convert to a GRange object
centro_telomere <- makeGRangesFromDataFrame(centro_telomere,
                                            keep.extra.columns=FALSE,
                                            ignore.strand=TRUE,
                                            seqinfo=NULL,
                                            seqnames.field="chr",
                                            start.field="start",
                                            end.field="end",
                                            starts.in.df.are.0based=FALSE)

# TCGA projects
TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

# remove cancer type which n < 100
TCGA_projects <- TCGA_projects[!(TCGA_projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", "THYM", "UCS", "UVM"))]


###################################################################################################
### Logistic regression to test whether age associates with increased/decreased possibility of focal regions to be gained or lost

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}


# remove centromere/telomere regions
remove_centro_telo <- function(regions){
  tmp <- as.character(regions$Peak.Limits)
  tmp2 <- c()
  for(i in tmp){
    tmp2 <- c(tmp2, unlist(strsplit(i, split = "[(]"))[1])
  }
  tmp_list <- list()
  for(i in 1:length(tmp2)){
    tmp_list[[i]] <- unlist(strsplit(tmp2[i], split = "[:,-]"))
  }
  tmp_df <- as.data.frame(do.call(rbind, tmp_list))
  colnames(tmp_df) <- c("chr", "start", "end")
  rownames(tmp_df) <- regions$Unique.Name
  
  tmp_df <- makeGRangesFromDataFrame(tmp_df,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)
  
  head(tmp_df)
  head(centro_telomere)
  
  tmp_df <- tmp_df[!(overlapsAny(tmp_df, centro_telomere)),]  # filter regions overlapping with centromeres or telomeres
  tmp_df <- as.data.frame(tmp_df)
  return(rownames(tmp_df))
}


# test association between age and focal-level CNAs for each region
test_asso <- function(name, project, df, match_df, GainOrLoss){
  region <- as.character(match_df[match_df$Unique.Name == name,]$Descriptor)
  region <- gsub(" ", "", region, fixed = TRUE)
  
  name_tmp <- unlist(strsplit(name, split = "Peak "))[2]
  peak_num <- unlist(strsplit(name_tmp, split = " - CN"))[1]  # number of Amp peaks or loss peaks we are working on (will use later for naming the result file)
  peak_num <- gsub(pattern = " ", replacement = "", peak_num, fixed = TRUE)
  
  df_row <- as.data.frame(df[rownames(df) == name,])
  df_row$patient <- rownames(df_row)
  rownames(df_row) <- NULL
  df_row <- df_row[,c(2,1)]
  colnames(df_row) <- c("patient", "value")
  
  # merge with age data
  tmp_df <- merge(df_row, clinical, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- glm(value ~ age , data = tmp_df, family = "binomial")
  summary(logit_fit)
  p_value <- formatC(as.numeric(summary(logit_fit)$coefficients[,4][2]), format = "e", digits = 2)
  coeff <- summary(logit_fit)$coefficients[,1][2]
  Z <- summary(logit_fit)$coefficients[,3][2]
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  
  write.csv(result, paste0("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Univariate_age_recurrent_focal_", GainOrLoss, 
                           "/", project, "_univariate_age_recurrent_", GainOrLoss, "_peak", peak_num, "_", region, ".csv", collapse = ""),  
            row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- name
  result_df$cancer_type <- project
  result_df$region <- region
  colnames(result_df) <- c("Unique.Name", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "cancer_type", "region")
  result_df <- result_df[,c("cancer_type", "Unique.Name", "region", "estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")]
  return(result_df)
}


# function to test gain/loss for each cancer type
age_focal <- function(project){
  
  print(paste0("working on : ", project))
  ### read significant results
  Sig_df <- read.delim(paste0("Analysis_results/CNAs/1_GISTIC2_CNAs/", project, "_GISTIC_mkCNV/all_lesions.conf_95.txt", collapse = ""))
  Sig_df$X <- NULL   # remove last column which has nothing but NAs
  Sig_df <- Sig_df[Sig_df$Amplitude.Threshold == "Actual Copy Change Given", ]   # keep only rows in section 2 (which has actual copy number changes)
  
  Sig_df <- Sig_df[Sig_df$Unique.Name %in% remove_centro_telo(Sig_df),]   # remove peaks overlapping with centromere or telomere
  rownames(Sig_df) <- Sig_df$Unique.Name
  
  # matching df for peak name and regions
  match_df <- Sig_df[,c("Unique.Name", "Descriptor")]
  rownames(match_df) <- NULL                 
  match_df$Descriptor <- unlist(lapply(as.character(match_df$Descriptor), gsub, pattern = " ", replacement = "", fixed = TRUE))
  
  df <- Sig_df[, -(7:9)]    # remove some columns
  df <- df[,-(3:6)]   # remove some columns to keep only columns for samples
  
  df_gain <- df[grepl('Amp', df[,1]),]    # gain rows
  n_gain <- nrow(df_gain)   # number of gain peaks
  
  df_loss <- df[grepl('Del', df[,1]),]    # loss rows
  n_loss <- nrow(df_loss)     # number of loss peaks
  
  # check number of gain samples
  num_gain <- c()
  for(i in 1:nrow(df_gain)){
    num_gain <- c(num_gain, sum(as.numeric(df_gain[i,3:ncol(df_gain)]) > 0.25))
  }
  
  # check number of lost samples
  num_loss <- c()
  for(i in 1:nrow(df_loss)){
    num_loss <- c(num_loss, sum(as.numeric(df_loss[i,3:ncol(df_loss)]) < -0.25))
  }
  
  # freq
  freq_gain <- num_gain/length(3:ncol(df_gain))
  freq_loss <- num_loss/length(3:ncol(df_loss))
  
  ### keep only focal regions that were altered more than 5% of the samples
  df_gain <- df_gain[freq_gain > 0.05, ]
  df_loss <- df_loss[freq_loss > 0.05, ]
  
  df_gain$Unique.Name <- NULL    # remove the first column (Unique.Name)
  df_loss$Unique.Name <- NULL    # remove the first column (Unique.Name)
  
  df_gain$Descriptor <- unlist(lapply(as.character(df_gain$Descriptor), gsub, pattern = " ", replacement = "", fixed = TRUE))
  df_loss$Descriptor <- unlist(lapply(as.character(df_loss$Descriptor), gsub, pattern = " ", replacement = "", fixed = TRUE))
  
  gain_regions <- df_gain$Descriptor
  loss_regions <- df_loss$Descriptor
  
  df_gain$Descriptor <- NULL
  df_loss$Descriptor <- NULL
  
  colnames(df_gain) <- unlist(lapply(colnames(df_gain), clean_id))
  colnames(df_loss) <- unlist(lapply(colnames(df_loss), clean_id))
  
  ### convert dataframe to TRUE, FALSE
  df_gain <- ifelse(df_gain > 0.25, TRUE, FALSE)   # check if gain or not
  df_loss <- ifelse(df_loss < -0.25, TRUE, FALSE)   # check if loss or not
  
  #### test association
  if(n_gain > 0){
    tmp_gain <- lapply(rownames(df_gain), test_asso, project = project, df = df_gain, match_df = match_df, GainOrLoss = "gain")
    tmp_gain <- do.call(rbind, tmp_gain)
    
    # calculate q values
    tmp_gain$p.value <- as.numeric(as.character(tmp_gain$p.value))
    tmp_gain$q.value <- p.adjust(tmp_gain$p.value, method = "BH")
    tmp_gain$GainOrLoss <- rep("gain", nrow(tmp_gain))
  }
  
  if(n_loss > 0){
    tmp_loss <- lapply(rownames(df_loss), test_asso, project = project, df = df_loss, match_df = match_df, GainOrLoss = "loss")
    tmp_loss <- do.call(rbind, tmp_loss)
    
    # calculate q values
    tmp_loss$p.value <- as.numeric(as.character(tmp_loss$p.value))
    tmp_loss$q.value <- p.adjust(tmp_loss$p.value, method = "BH")
    tmp_loss$GainOrLoss <- rep("loss", nrow(tmp_loss))
  }
  
  if(n_gain > 0 & n_loss > 0){
    tmp <- rbind(tmp_gain, tmp_loss)
    print(paste0("finish : ", project))
    return(tmp)
  } else if (n_gain > 0){
    return(tmp_gain)
  } else if (n_loss > 0){
    return(tmp_loss)
  }
}

### call function #################################################################################
all_cancers_results <- lapply(TCGA_projects, age_focal)
all_cancers_results <- do.call(rbind, all_cancers_results)
all_cancers_results$Threshold <- ifelse(all_cancers_results$q.value < 0.05, TRUE, FALSE)   # check if gain or not

rownames(all_cancers_results) <- NULL
head(all_cancers_results)

all_cancers_results$odds <- exp(all_cancers_results$estimate)
all_cancers_results$odds_conf.low <- exp(all_cancers_results$conf.low)
all_cancers_results$odds_conf.high <- exp(all_cancers_results$conf.high)
all_cancers_results <- all_cancers_results[,c("cancer_type", "Unique.Name", "region", "estimate", "std.error",
                                              "conf.low", "conf.high", "statistic", "odds", "odds_conf.low", "odds_conf.high",
                                              "p.value", "q.value", "GainOrLoss", "Threshold")]
write.csv(all_cancers_results, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Univariate_age_recurrent_focal.csv", row.names = FALSE)

sig_result <- all_cancers_results[all_cancers_results$Threshold == TRUE,]
write.csv(sig_result, "Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Univariate_age_recurrent_focal_sig.csv", row.names = FALSE)
