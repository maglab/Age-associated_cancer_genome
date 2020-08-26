# check missing clinical data and save useful clinical data in new files
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(stringr)

# clinical file from all_clin_indexed.csv
clin <- read.csv("Data/all_clin_indexed.csv")

# samples included in the study
samples <- read.csv("Data/samples_in_ASCAT_and_seg.csv")
samples <- as.character(samples$patient)

TCGA_projects <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                   "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                   "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

### ACC ###########################################################################################
project <- "ACC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples
  
age.na <- sum(is.na(clinical$primary_pathology_age_at_initial_pathologic_diagnosis))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))

clinical <- clinical[,c("bcr_patient_barcode", "primary_pathology_age_at_initial_pathologic_diagnosis", "gender", "race_list", "stage_event_pathologic_stage")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)

### BLCA ##########################################################################################
project <- "BLCA"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))
subtype.na <- sum(is.na(clinical$diagnosis_subtype))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade", 
                         "diagnosis_subtype", "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade", 
                        "subtype", "smoking_history")

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### BRCA ##########################################################################################
project <- "BRCA"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
ER.na <- sum(is.na(clinical$breast_carcinoma_estrogen_receptor_status))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "breast_carcinoma_estrogen_receptor_status")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "ER_status")

pathologic <- clinical$pathologic_stage
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$pathologic_stage <- as.character(pathologic)

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### CESC ##########################################################################################
project <- "CESC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
figo.stage.na <- sum(is.na(clinical$stage_event_clinical_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_clinical_stage", "neoplasm_histologic_grade", 
                         "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "figo_stage", 
                        "histologic_grade", "smoking_history")

figo <- as.character(clinical$figo_stage)
figo[figo %in% c("Stage I", "Stage IA", "Stage IA2", "Stage IB", "Stage IB1", "Stage IB2")] <- "Stage I"
figo[figo %in% c("Stage II", "Stage IIA", "Stage IIA1", "Stage IIA2", "Stage IIB")] <- "Stage II"
figo[figo %in% c("Stage III", "Stage IIIA", "Stage IIIB")] <- "Stage III"
figo[figo %in% c("Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$figo_stage <- figo

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### CHOL ##########################################################################################
project <- "CHOL"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))

risk_factor <- as.character(clinical$history_hepato_carcinoma_risk_factors)
clinical$smoking_history <- unlist(lapply(risk_factor, str_detect, pattern = "Smoking"))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "smoking_history")]

colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "smoking_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)

### COAD ##########################################################################################
project <- "COAD"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
subtype.na <- sum(is.na(clinical$histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "histological_type")]

colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "subtype")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### DLBC ##########################################################################################
project <- "DLBC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
subtype.na <- sum(is.na(clinical$histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", "histological_type")]

colnames(clinical) <- c("patient", "age", "gender", "race", "subtype")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)

### ESCA ##########################################################################################
project <- "ESCA"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$primary_pathology_neoplasm_histologic_grade))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))
alcohol.na <- sum(is.na(clinical$alcohol_history_documented))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "primary_pathology_neoplasm_histologic_grade", 
                         "tobacco_smoking_history", "alcohol_history_documented")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade", 
                        "smoking_history", "alcohol_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### GBM ##########################################################################################
project <- "GBM"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list")]

colnames(clinical) <- c("patient", "age", "gender", "race")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### HNSC ##########################################################################################
project <- "HNSC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))
alcohol.na <- sum(is.na(clinical$alcohol_history_documented))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade", 
                         "tobacco_smoking_history", "alcohol_history_documented")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade", 
                        "smoking_history", "alcohol_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage IVA", "Stage IVB", "Stage IVC")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "G3"] <- "G3.G4"
clinical$histologic_grade[clinical$histologic_grade == "G4"] <- "G3.G4"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### KICH ##########################################################################################
project <- "KICH"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "smoking_history")

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### KIRC ##########################################################################################
project <- "KIRC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade", 
                         "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade", 
                        "smoking_history")

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "G3"] <- "G3.G4"
clinical$histologic_grade[clinical$histologic_grade == "G4"] <- "G3.G4"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### KIRP ##########################################################################################
project <- "KIRP"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "smoking_history")

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### LAML ##########################################################################################
project <- "LAML"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list")]
colnames(clinical) <- c("patient", "age", "gender", "race")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### LGG ##########################################################################################
project <- "LGG"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", "neoplasm_histologic_grade")]
colnames(clinical) <- c("patient", "age", "gender", "race", "histologic_grade")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### LIHC ##########################################################################################
project <- "LIHC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))

risk_factor <- as.character(clinical$history_hepato_carcinoma_risk_factors)
clinical$alcohol_history <- unlist(lapply(risk_factor, str_detect, pattern = "Alcohol consumption"))
clinical$Hepatitis <- unlist(lapply(risk_factor, str_detect, pattern = "Hepatitis"))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade", 
                         "alcohol_history", "Hepatitis")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade", 
                        "alcohol_history", "Hepatitis")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "G3"] <- "G3.G4"
clinical$histologic_grade[clinical$histologic_grade == "G4"] <- "G3.G4"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### LUAD ##########################################################################################
project <- "LUAD"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "smoking_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### LUSC ##########################################################################################
project <- "LUSC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))

smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "tobacco_smoking_history")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "smoking_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### MESO ##########################################################################################
project <- "MESO"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
subtype.na <- sum(is.na(clinical$primary_pathology_histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "primary_pathology_histological_type")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "subtype")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### OV ##########################################################################################
project <- "OV"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
figo.stage.na <- sum(is.na(clinical$stage_event_clinical_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_clinical_stage", "neoplasm_histologic_grade")]
colnames(clinical) <- c("patient", "age", "gender", "race", "figo_stage", "histologic_grade")

figo <- as.character(clinical$figo_stage)
figo[figo %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC")] <- "Stage I"
figo[figo %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
figo[figo %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$figo_stage <- figo

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "G3"] <- "G3.G4"
clinical$histologic_grade[clinical$histologic_grade == "G4"] <- "G3.G4"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### PAAD ##########################################################################################
project <- "PAAD"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))
smoking_history.na <- sum(is.na(clinical$tobacco_smoking_history))
alcohol.na <- sum(is.na(clinical$alcohol_history_documented))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade", 
                         "tobacco_smoking_history", "alcohol_history_documented")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", 
                        "histologic_grade", "smoking_history", "alcohol_history")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"

clinical$pathologic_stage <- pathologic

clinical$smoking_history[clinical$smoking_history == 1] <- "Lifelong Non-smoker"
clinical$smoking_history[clinical$smoking_history == 2] <- "Current smoker"
clinical$smoking_history[clinical$smoking_history == 3] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 4] <- "Current reformed smoker"
clinical$smoking_history[clinical$smoking_history == 5] <- "Current reformed smoker"

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "G3"] <- "G3.G4"
clinical$histologic_grade[clinical$histologic_grade == "G4"] <- "G3.G4"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### PCPG ##########################################################################################
project <- "PCPG"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list")]
colnames(clinical) <- c("patient", "age", "gender", "race")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### PRAD ##########################################################################################
project <- "PRAD"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index", "primary_gleason_grade", "secondary_gleason_grade")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")
primary_gleason <- as.character(clinical$primary_gleason_grade)
primary_gleason <- as.numeric(unlist(lapply(primary_gleason, gsub, pattern = "Pattern ", replacement = "")))
clinical$primary_gleason_grade <- primary_gleason
secondary_gleason <- as.character(clinical$secondary_gleason_grade)
secondary_gleason <- as.numeric(unlist(lapply(secondary_gleason, gsub, pattern = "Pattern ", replacement = "")))
clinical$secondary_gleason_grade <- secondary_gleason
clinical$gleason_score <- clinical$primary_gleason_grade + clinical$secondary_gleason_grade

gleason_score <- c()
for(i in 1:nrow(clinical)){
  if(clinical$gleason_score[i] <= 6){
    gleason_score <- c(gleason_score, "<=6")
  } else if(clinical$gleason_score[i] >= 8){
    gleason_score <- c(gleason_score, ">=8")
  } else if(clinical$gleason_score[i] == 7 & clinical$primary_gleason_grade[i] == 3){
    gleason_score <- c(gleason_score, "7(3+4)")
  } else if(clinical$gleason_score[i] == 7 & clinical$primary_gleason_grade[i] == 4){
    gleason_score <- c(gleason_score, "7(4+3)")
  }
}
clinical$gleason_score <- gleason_score

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
gleason_score.na <- sum(is.na(clinical$gleason_score))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", "gleason_score")]
colnames(clinical) <- c("patient", "age", "gender", "race", "gleason_score")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### READ ##########################################################################################
project <- "READ"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
subtype.na <- sum(is.na(clinical$histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "histological_type")]

colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "subtype")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
pathologic[pathologic %in% c("Stage IV", "Stage IVA")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### SARC ##########################################################################################
project <- "SARC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
subtype.na <- sum(is.na(clinical$primary_pathology_histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index", "gender", "race_list", "primary_pathology_histological_type")]

colnames(clinical) <- c("patient", "age", "gender", "race", "subtype")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### SKCM ##########################################################################################
project <- "SKCM"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
clark.na <- sum(is.na(clinical$melanoma_clark_level_value))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index", "gender", "race_list", 
                         "stage_event_pathologic_stage", "melanoma_clark_level_value")]

colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "clark_value")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### STAD ##########################################################################################
project <- "STAD"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "neoplasm_histologic_grade")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "histologic_grade")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### TGCT ##########################################################################################
project <- "TGCT"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))

histology <- as.character(clinical$primary_pathology_histology_list)
subtype <- c()
for(i in histology){
  if(str_detect(i, pattern = "Non-Seminoma")){
    subtype <- c(subtype, "Non-Seminoma")
  } else if(str_detect(i, pattern = "Seminoma")){
    subtype <- c(subtype, "Seminoma")
  }
}

clinical$subtype <- subtype

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "subtype")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "subtype")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### THCA ##########################################################################################
project <- "THCA"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))
subtype.na <- sum(is.na(clinical$histological_type))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_pathologic_stage", "histological_type")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage", "subtype")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC")] <- "Stage IV"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### THYM ##########################################################################################
project <- "THYM"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list")]
colnames(clinical) <- c("patient", "age", "gender", "race")

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### UCEC ##########################################################################################
project <- "UCEC"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
figo.stage.na <- sum(is.na(clinical$stage_event_clinical_stage))
grade.na <- sum(is.na(clinical$neoplasm_histologic_grade))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_clinical_stage", "neoplasm_histologic_grade")]
colnames(clinical) <- c("patient", "age", "gender", "race", "figo_stage", "histologic_grade")

figo <- as.character(clinical$figo_stage)
figo[figo %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC")] <- "Stage I"
figo[figo %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
figo[figo %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IIIC1", "Stage IIIC2")] <- "Stage III"
figo[figo %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$figo_stage <- figo

clinical$histologic_grade <- as.character(clinical$histologic_grade)
clinical$histologic_grade[clinical$histologic_grade == "High Grade"] <- "G3"

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### UCS ##########################################################################################
project <- "UCS"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
figo.stage.na <- sum(is.na(clinical$stage_event_clinical_stage))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", 
                         "stage_event_clinical_stage")]
colnames(clinical) <- c("patient", "age", "gender", "race", "figo_stage")

figo <- as.character(clinical$figo_stage)
figo[figo %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC")] <- "Stage I"
figo[figo %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
figo[figo %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IIIC1", "Stage IIIC2")] <- "Stage III"
figo[figo %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"

clinical$figo_stage <- figo

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)


### UVM ###########################################################################################
project <- "UVM"
clinical <- read.csv(paste0("Data/clinical_XML/TCGA-", project, "_clinical_XML.csv", collapse = ""), na.strings=c("","NA"))
clinical <- clinical[clinical$bcr_patient_barcode %in% samples,]
clinical <- clinical[!duplicated(clinical), ]   # remove duplicated rows

num_samples <- nrow(clinical)   # numbers of samples

# get age from all_clin_index
clin_tmp <- clin[,c("submitter_id", "age_at_index")]
clinical <- merge(clinical, clin_tmp, by.x = "bcr_patient_barcode", by.y = "submitter_id")

age.na <- sum(is.na(clinical$age_at_index))
gender.na <- sum(is.na(clinical$gender))
race.na <- sum(is.na(clinical$race_list))
stage.na <- sum(is.na(clinical$stage_event_pathologic_stage))

clinical <- clinical[, c("bcr_patient_barcode", "age_at_index","gender", "race_list", "stage_event_pathologic_stage")]
colnames(clinical) <- c("patient", "age", "gender", "race", "pathologic_stage")

pathologic <- as.character(clinical$pathologic_stage)
pathologic[pathologic %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
pathologic[pathologic %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"

clinical$pathologic_stage <- pathologic

write.csv(clinical, paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""), row.names = FALSE)

###################################################################################################
######################### Number of samples for each group of variables ###########################
num_samples <- function(project){
  clinical <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  clinical <- clinical[,-c(1,2)]
  print(apply(clinical, 2, table, exclude = NULL))
}


samples <- lapply(TCGA_projects,num_samples)
names(samples) <- TCGA_projects
samples

### get only age, sex, race (will be used in PANCAN analysis)
get_AgeSexRace <- function(project){
  clinical <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  clinical$cancer_type <- rep(project, nrow(clinical))
  df <- clinical[,c("patient", "cancer_type", "age", "gender", "race")]
  return(df)
}

clinical_df <- lapply(TCGA_projects, get_AgeSexRace)

result_df <- do.call(rbind, clinical_df)

write.csv(result_df, "Data/all_clin_XML.csv", row.names = FALSE)


