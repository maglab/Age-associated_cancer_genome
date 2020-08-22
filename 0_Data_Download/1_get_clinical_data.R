# Download clinical data
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/Data")

# Get clinical data using TCGAbiolinks
library(TCGAbiolinks)

# This code will get all clinical indexed data from TCGA
library(data.table)
library(dplyr)
library(regexPipes)

clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist(fill = TRUE)
readr::write_csv(clinical,path = paste0("Data/", "all_clin_indexed.csv", collapse = ""))

# get clinical XML
TCGA_projects <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% sort

for(project in TCGA_projects){
  query <- GDCquery(project = project, 
                    data.category = "Clinical", 
                    file.type = "xml")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  write.csv(clinical, paste0("clinical_XML/", project, "_clinical_XML.csv", collapse = ""), row.names = FALSE)
  
}

# MSI status
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/COAD_MSI_Status.csv", row.names = FALSE)

query <- GDCquery(project = "TCGA-READ", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/READ_MSI_Status.csv", row.names = FALSE)

query <- GDCquery(project = "TCGA-ESCA", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/ESCA_MSI_Status.csv", row.names = FALSE)

query <- GDCquery(project = "TCGA-STAD", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/STAD_MSI_Status.csv", row.names = FALSE)

query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/PAAD_MSI_Status.csv", row.names = FALSE)

query <- GDCquery(project = "TCGA-UCEC", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test")
GDCdownload(query)
msi_results <- GDCprepare_clinic(query, "msi")
msi_results %>% DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
write.csv(msi_results, "Data/MSI_Status/UCEC_MSI_Status.csv", row.names = FALSE)
