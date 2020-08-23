### download gene expression data
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/Data/TCGA_Gene_Expression")

library(TCGAbiolinks)
library(SummarizedExperiment)

TCGA_projects <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA",
                   "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC",
                   "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
                   "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")

download_expression <- function(project){
  cancer_type <- unlist(strsplit(project, split = "-"))[2]
  query <- GDCquery(project = project, 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query)
  data <- GDCprepare(query)
  
  write.csv(assay(data), paste0(project, "_normalised_expression.csv", collapse = ""))
}

lapply(TCGA_projects, download_expression)


### Download expression from Illumina GA platform
download_expression_GA <- function(project){
  cancer_type <- unlist(strsplit(project, split = "-"))[2]
  query <- GDCquery(project = project, 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina GA", 
                    file.type = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query)
  data <- GDCprepare(query)
  
  write.csv(assay(data), paste0(project, "_normalised_expression_GA.csv", collapse = ""))
}

lapply(c("TCGA-COAD", "TCGA-READ", "TCGA-UCEC"), download_expression_GA)


