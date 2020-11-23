### POLE and POLD1 multiple logistic regression analysis with age
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(logistf)
library(ggplot2)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

projects <- projects[!(projects %in% c("DLBC", "PCPG", "TGCT", "THYM", "UVM"))] # remove projects with no mut in POLE or POLD1

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)


# sub-function to select a logistic regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- mut ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- mut ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- mut ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- mut ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- mut ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- mut ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- mut ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- mut ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- mut ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- mut ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- mut ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- mut ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- mut ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- mut ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- mut ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- mut ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- mut ~ age + purity + gender + pathologic_stage
  }
  return(model)
}

# function to test association between age and mutation
test_age_mut_gene <- function(project, gene, purity){
  
  print(paste("Working on: ", project, ";", gene))
  
  ### clinical data
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  head(clin)
  
  ### read gene-sample mutation table
  mut_df <- read.csv(paste0("Analysis_results/Mutations/6_POLE_POLD1/", project, "_pol_mutations.csv", collapse = ""))
  row.names(mut_df) <- mut_df$X
  mut_df$X <- NULL
  
  ### select model
  model <- model_selection(project)
  
  # get value for gene of interest
  if(gene == "POLE" & !(gene %in% colnames(mut_df))){
    return(NULL)
  } else if(gene == "POLD1" & !(gene %in% colnames(mut_df))){
    return(NULL)
  }
  
  values_gene <- mut_df[gene]
  
  if(sum(values_gene[1]) > 0.05*nrow(values_gene)){ # work on projects that contain POLE or POLD1 mutations in > 5% of samples
    values_gene$patient <- rownames(values_gene)
    rownames(values_gene) <- NULL
    colnames(values_gene) <- c("mut", "patient")
    
    # merge with age data
    df <- merge(values_gene, clin, by.x = "patient", by.y = "patient")
    
    # merge with purity
    df <- merge(df, purity, by.x = "patient", by.y = "patient")
    
    # write source data
    if(project == "UCEC"){
      write.csv(df, paste0("Source_Data/Fig_4d_", gene, ".csv", collapse = ""), row.names = FALSE)
    }else{
      write.csv(df, paste0("Source_Data/Supplementary_Fig_8b_", project, "_", gene, ".csv", collapse = ""), row.names = FALSE)
    }
    
    # logistic regression
    logit_fit <- logistf(formula = model, data = df, family = "binomial")
    
    summary(logit_fit)
    
    result <- broomExtra::tidy_parameters(logit_fit)
    
    write.csv(result, paste0("Analysis_results/Mutations/6.1_POLE_POLD1_results/", project, "_multivariate_age_mut_", gene, ".csv", collapse = ""), row.names = FALSE)
    
    result_df <- as.data.frame(result[result$term == "age",])
    result_df$term <- gene
    result_df$cancer_type <- project
    colnames(result_df) <- c("gene", "estimate", "std.error", "conf.low", "conf.high", "statistics", "df.error", "p.value", "cancer_type")
    result_df <- result_df[,c("cancer_type", "gene", "estimate", "std.error", "conf.low", "conf.high", "statistics", "df.error", "p.value")]
    
    ### plot
    df$mut_1 <- ifelse(df$mut == 1, TRUE, FALSE)
    
    ### Fig. 4d and Supplementary Fig. 8b
    my_label <- paste0("p = ", round(result_df$p.value, 4))
    pdf(paste0("Analysis_results/Mutations/6.1_POLE_POLD1_results/", project, "_age_", gene,"_mut.pdf", collapse = ""),
        width = 3, height = 4, useDingbats=FALSE) 
    p <- ggplot(aes(x = mut_1, y = age, fill = mut_1), data = df) + 
      geom_violin(trim = FALSE, scale = "width") + 
      geom_boxplot(width = 0.4, fill = "white") +
      ggtitle(paste0(project, ": ", gene, " mutation", collapse = "")) +
      xlab("mutation") +
      ylab("Age at diagnosis") +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size=14,face="bold"),
            axis.title.y = element_text(size=14,face="bold"),
            legend.title = element_blank(),
            legend.position = "none",
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      annotate("label", x=-Inf, y = Inf, size = 5,
               label = my_label, hjust=0, vjust=1)
    print(p)
    dev.off()
    return(result_df)
  } else {
    return(NULL)
  }
}

POLE_results <- lapply(projects, test_age_mut_gene, gene = "POLE", purity = purity)
POLE_results[sapply(POLE_results, is.null)] <- NULL   # remove NULL
POLE_results <- do.call(rbind, POLE_results)
write.csv(POLE_results, "Analysis_results/Mutations/6.1_POLE_POLD1_results/POLE_age_mut.csv")

POLD1_results <- lapply(projects, test_age_mut_gene, gene = "POLD1", purity = purity)
POLD1_results[sapply(POLD1_results, is.null)] <- NULL   # remove NULL
POLD1_results <- do.call(rbind, POLD1_results)
write.csv(POLD1_results, "Analysis_results/Mutations/6.1_POLE_POLD1_results/POLD1_age_mut.csv")

