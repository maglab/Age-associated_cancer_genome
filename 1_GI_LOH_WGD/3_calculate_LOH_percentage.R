### Script for calculating percent genomic LOH ####################################################
# input files: ASCAT profiles
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

### Chromosome region file ########################################################################
chrom <- read.table("Data/hg19_chrom_length.txt", header = TRUE)  # from UCSC genome browser

# calculate genome length excluding X and Y
chrom_no_XY <- chrom[!(chrom$chr %in% c("chrX", "chrY")),]
genome_length <- 0
for(i in 1:nrow(chrom_no_XY)){
  genome_length <- genome_length + chrom_no_XY[i,]$length_p + chrom_no_XY[i,]$length_q
}

genome_length

### helper functions ##############################################################################
# 1) add arm to regions
add_arm <- function(region){
   chr <- region$chr  # chromosome
   chr <- paste0("chr", chr, collapse = "")
   chrom_chr <- chrom[chrom$chr == chr,]   # get chromosome info from chrom dataframe
   if(region$startpos < chrom_chr$p_end & region$endpos < chrom_chr$p_end){return("p")}
   else if(region$startpos > chrom_chr$p_end & region$endpos < chrom_chr$q_end){return("q")}
   else {return("both")}
}

# 2) determine if it contain LOH
# if the regions is LOH, return region length
# if not, return FALSE
check_LOH <- function(region, arm){
  if(region$nMajor == 0 & region$nMinor == 0){
    return(FALSE)
  } else if(region$nMajor != 0 & region$nMinor != 0){
    return(FALSE)
  } else if(region$nMajor == 0 | region$nMinor == 0){
    region_length <- region$endpos - region$startpos  # length of affected region
    chr <- paste0("chr", region$chr, collapse = "")   # chrom
    
    # get chromosome length
    chromlength <- chrom[chrom$chr == chr,]
    if(arm == "p"){
      chromlength <- chromlength$length_p
    } else if(arm == "p"){
      chromlength <- chromlength$length_q
    } else{
      chromlength <- min(c(chromlength$length_p, chromlength$length_q))
    }
    
  }
    return(region_length)
}

### main function #################################################################################

percent_LOH <- function(ASCAT_file, genome_length){
  
  ASCAT_file <- read.table(ASCAT_file, header = TRUE)
  ASCAT_file <- ASCAT_file[!(ASCAT_file$chr %in% c("X", "Y")),]   # remove X and Y
  
  all_LOH_regions <- 0      # for LOH length
  
  for(i in 1:nrow(ASCAT_file)){
    region <- ASCAT_file[i,]
    arm <- add_arm(region)
    LOH_status <- check_LOH(region, arm)
    if(is.integer(LOH_status)){
      all_LOH_regions <- all_LOH_regions + LOH_status
    }
  }
  
  return((all_LOH_regions/genome_length) * 100)
}


### work on all ASCAT files #######################################################################
# list files
folders <-list.dirs("Data/ASCAT_TCGA_filtered")
folders <- folders[folders != "Data/ASCAT_TCGA_filtered"]

result_list <- list()

for(folder in folders){
  project <- unlist(strsplit(folder, split = "\\/"))[3]
  print(paste0("working on: ", project))
  
  files <- list.files(folder, pattern = ".txt", full.names=TRUE)
  for(f in files){
    file_name <- unlist(strsplit(f, split = "\\/"))[4]
    file_name <- unlist(strsplit(file_name, split = "[.]"))[1]
    LOH <- percent_LOH(f, genome_length = genome_length)
    result_list[[file_name]] <- c(project, LOH)
  }
}

df <- do.call(rbind, result_list)
df <- as.data.frame(df)
df$patient <- rownames(df)
rownames(df) <- NULL

colnames(df) <- c("cancer_type", "percent_LOH", "patient")
df <- df[,c(3,1,2)]

# clean id
clean_id <- function(patient){
  tmp <- unlist(strsplit(patient, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

df$patient <- unlist(lapply(as.character(df$patient), clean_id))

write.table(df, "Data/TCGA.LOH_percentage.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)


