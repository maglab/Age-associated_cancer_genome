### script to download somatic copy number file from TCGA using TCGAbiolinks
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/Data/TCGA_segfiles/")

library(TCGAbiolinks)
library(dplyr)

projects <- getGDCprojects()$project_id %>% regexPipes::grep("TCGA",value=T) %>% sort

get_seg <- function(project){
  query <- GDCquery(project = project, 
                    data.category = "Copy number variation",
                    file.type = "nocnv_hg19.seg",
                    legacy = TRUE)
  GDCdownload(query)
  
  # check duplicated samples
  # remove duplicated samples by keeping the first one
  if(sum(duplicated(query$results[[1]]$cases)) > 0){
    print(paste0("Duplicated samples in: ", project, " ; Number of duplicated samples: ", sum(duplicated(query$results[[1]]$cases)), collapse = ""))
    query$results[[1]] <- query$results[[1]][!duplicated(query$results[[1]]$cases),]
  }
  data <- GDCprepare(query)
  write.table(data, paste0(project, "_noCNV.txt", collapse = ""), quote = FALSE, row.names = FALSE, sep = "\t")
  
  return(project)
}

lapply(projects, get_seg)
