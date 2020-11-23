### 6 substitution classes with age for each cancer

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(maftools)
library(ggplot2)
library(broom)
library(reshape2)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

########## Plot mut fraction contribution with age for each cancer type (simple linear regression) ##############################
### linear regression for each mutation class
mut_age_class_fraction_lm <- function(class, mut_fraction_df, clinical, project){
  
  print(paste0("Working on: ", class))
  
  if(class == "C>A"){
    my_col <- "#F8766D"
  } else if(class == "C>G"){
    my_col <- "#B79F00"
  } else if(class == "C>T"){
    my_col <- "#00BFC4"
  } else if(class == "T>C"){
    my_col <- "#00BA38"
  } else if(class == "T>A"){
    my_col <- "#619CFF"
  } else if(class == "T>G"){
    my_col <- "#F564E3"
  }
  
  tmp_df <- mut_fraction_df[,c("Tumor_Sample_Barcode", class)]
  colnames(tmp_df) <- c("patient", "fraction")
  
  tmp_df <- merge(tmp_df, clinical, by.x = "patient", by.y = "patient")
  
  # linear regression age and mutational burden
  lm_fit <- lm(fraction ~ age, data=tmp_df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared,2)
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  pdf(paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/", project, "_univariate_fraction_",class,".pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE)
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  p <- ggplot(data = tmp_df, aes(x = age, y = fraction)) + 
    geom_point(color = my_col) +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle(paste0(project, ": ", class)) +
    xlab("Age at diagnosis") +
    ylab("Fraction contribution") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face="bold"),
          panel.background = element_blank(),
          #panel.border = element_rect(linetype = "solid", fill = NA),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 6,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()

  result_df <- tidy(lm_fit)
  
  write.csv(result_df, paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/",project, "_univariate_fraction_",class,".csv", collapse = ""), row.names = FALSE)
  
  out <- result_df[result_df$term == "age",]
  out$term <- class
  colnames(out) <- c("class", "estimate", "std.error", "statistic", "p.value")
  return(out)
}


### function to test lm for each project univariate
plot_mut_fraction_age <- function(project, clinical){
  
  print(paste0("Working on: ", project))
  
  # read maf
  maf <- read.maf(maf = paste0("Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""))
  
  maf@variant.classification.summary
  titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
  mut_fraction_df <- as.data.frame(titv$fraction.contribution)
  mut_fraction_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_fraction_df$Tumor_Sample_Barcode), clean_id))
  
  classes <- c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")
  
  ### test and plot
  result_project <- lapply(classes, mut_age_class_fraction_lm, mut_fraction_df, clinical, project = project)
  result_project <- do.call(rbind, result_project)
  result_project <- as.data.frame(result_project)
  result_project$q.value <- p.adjust(result_project$p.value, method = "BH")
  result_project$Sig <- ifelse(result_project$q.value < 0.05, TRUE, FALSE)
  write.csv(result_project, paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/",project, "_univariate_fraction_all_classes.csv", collapse = ""), row.names = FALSE)
  
  return(paste0("Finish: ", project, ", ", nrow(mut_fraction_df)))
}

lapply(projects, plot_mut_fraction_age, clinical = clinical)


#### Check sig all projects 
read_project <- function(project){
  tmp <- read.csv(paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/",project, "_univariate_fraction_all_classes.csv", collapse = ""))
  tmp$project <- rep(project, nrow(tmp))
  tmp <- tmp[,c("project", "class", "estimate", "std.error", "statistic", "p.value", "q.value", "Sig")]
  return(tmp)
}

df <- lapply(projects, read_project)
df <- do.call(rbind, df)
df$direction <- ifelse(df$estimate > 0, "up", "down")
write.csv(df, "Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_univariate_fraction_all.csv", row.names = FALSE)

sig_df <- df[df$Sig == TRUE,]

write.csv(sig_df, "Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_sig_univariate_fraction_all.csv", row.names = FALSE)


########## Multiple linear regression ###########################################################################################
### read purity
purity <- read.table("Data/TCGA.purity.txt", header = TRUE)

### read sig. result from simple regression
sig_uni <- read.csv("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_sig_univariate_fraction_all.csv")
projects <- as.character(unique(sig_uni$project))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

# sub-function to select a regression model depends on cancer type
model_selection <- function(project){
  if(project %in% c("ACC")){
    model <- fraction ~ age + purity + gender + pathologic_stage
  } else if(project %in% c("BLCA")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + histologic_grade + subtype + smoking_history
  } else if(project %in% c("BRCA")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + ER_status
  } else if(project %in% c("CESC")){
    model <- fraction ~ age + purity + figo_stage + histologic_grade
  } else if(project %in% c("CHOL", "KICH", "LUAD")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + smoking_history
  } else if(project %in% c("COAD", "READ")){
    model <- fraction ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("DLBC", "SARC")){
    model <- fraction ~ age + purity + gender + race + subtype
  } else if(project %in% c("ESCA")){
    model <- fraction ~ age + purity + gender + histologic_grade + alcohol_history
  } else if(project %in% c("GBM", "LAML", "PCPG", "THYM")){
    model <- fraction ~ age + purity + gender + race
  } else if(project %in% c("HNSC")){
    model <- fraction ~ age + purity + gender + race + histologic_grade + smoking_history + alcohol_history
  } else if(project %in% c("KIRC")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + histologic_grade
  } else if(project %in% c("KIRP", "SKCM")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage
  } else if(project %in% c("LGG")){
    model <- fraction ~ age + purity + gender + race + histologic_grade
  } else if(project %in% c("LIHC")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history + Hepatitis
  } else if(project %in% c("LUSC")){
    model <- fraction ~ age + purity + gender + pathologic_stage + smoking_history
  } else if(project %in% c("MESO")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + subtype
  } else if(project %in% c("OV", "UCEC")){
    model <- fraction ~ age + purity + race + figo_stage + histologic_grade
  } else if(project %in% c("PAAD")){
    model <- fraction ~ age + purity + gender + race + pathologic_stage + histologic_grade + alcohol_history
  } else if(project %in% c("PRAD")){
    model <- fraction ~ age + purity + race + gleason_score
  } else if(project %in% c("STAD")){
    model <- fraction ~ age + purity + gender + pathologic_stage + histologic_grade
  } else if(project %in% c("TGCT")){
    model <- fraction ~ age + purity + race + pathologic_stage + subtype
  } else if(project %in% c("THCA")){
    model <- fraction ~ age + purity + gender + pathologic_stage + subtype
  } else if(project %in% c("UCS")){
    model <- fraction ~ age + purity + race + figo_stage
  } else if(project %in% c("UVM")){
    model <- fraction ~ age + purity + gender + pathologic_stage
  }
  return(model)
}


########## Plot mut fraction contribution with age for each cancer type (multiple regression) ###################################
### multiple linear regression for each mutation class
mut_age_class_fraction_lm <- function(class, mut_fraction_df, clin, purity, project, model){
  
  print(paste0("Working on: ", class))
  
  if(class == "C>A"){
    my_col <- "#F8766D"
  } else if(class == "C>G"){
    my_col <- "#B79F00"
  } else if(class == "C>T"){
    my_col <- "#00BFC4"
  } else if(class == "T>C"){
    my_col <- "#00BA38"
  } else if(class == "T>A"){
    my_col <- "#619CFF"
  } else if(class == "T>G"){
    my_col <- "#F564E3"
  }
  
  tmp_df <- mut_fraction_df[,c("Tumor_Sample_Barcode", class)]
  colnames(tmp_df) <- c("patient", "fraction")
  
  tmp_df <- merge(tmp_df, clin, by.x = "patient", by.y = "patient")
  tmp_df <- merge(tmp_df, purity, by.x = "patient", by.y = "patient")
  
  # linear regression age and mutational burden
  lm_fit <- lm(formula = model, data=tmp_df)
  p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
  p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
  r_squared <- round(summary(lm_fit)$r.squared,2)
  coeff <- summary(lm_fit)$coefficients[,1][2]
  
  pdf(paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/", project, "_multivariate_fraction_",class,".pdf", collapse = ""), width = 6, height = 4.5, useDingbats=FALSE)
  my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
  p <- ggplot(data = tmp_df, aes(x = age, y = fraction)) + 
    geom_point(color = my_col) +
    geom_smooth(method = "lm", se = TRUE) +
    ggtitle(paste0(project, ": ", class)) +
    xlab("Age at diagnosis") +
    ylab("Fraction contribution") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15,face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face="bold"),
          panel.background = element_blank(),
          #panel.border = element_rect(linetype = "solid", fill = NA),
          axis.line = element_line(colour = "black")) +
    annotate("label", x=-Inf, y = Inf, size = 6,
             label = my_label, hjust=0, vjust=1)
  print(p)
  dev.off()
  
  result_df <- tidy(lm_fit)
  
  write.csv(result_df, paste0("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/",project, "_multivariate_fraction_",class,".csv", collapse = ""), row.names = FALSE)
  
  out <- result_df[result_df$term == "age",]
  out$term <- class
  colnames(out) <- c("class", "estimate", "std.error", "statistic", "p.value")
  return(out)
}


### function to test lm for each project univariate
plot_mut_fraction_age <- function(project, sig_uni, purity){
  
  ### get information
  df_of_interest <- sig_uni[sig_uni$project == project,]
  project <- unique(as.character(df_of_interest$project))
  classes <- as.character(df_of_interest$class)
  
  print(paste0("Working on: ", project))
  
  # read maf
  maf <- read.maf(maf = paste0("Data/MC3_maf_TCGA_projects/", project, "_mc3.maf", collapse = ""))
  
  maf@variant.classification.summary
  titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
  mut_fraction_df <- as.data.frame(titv$fraction.contribution)
  mut_fraction_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_fraction_df$Tumor_Sample_Barcode), clean_id))
  
  ### clinical
  clin <- read.csv(paste0("Data/clinical_XML_interest/TCGA-", project, "_clinical_XML_interest.csv", collapse = ""))
  
  ### model selection
  model <- model_selection(project)
  
  ### test and plot
  result_project <- lapply(classes, mut_age_class_fraction_lm, mut_fraction_df, clin = clin, purity = purity, project = project, model = model)
  result_project <- do.call(rbind, result_project)
  result_project <- as.data.frame(result_project)
  result_project$project <- rep(project, nrow(result_project))
  result_project <- result_project[,c("project", "class", "estimate", "std.error", "statistic", "p.value")]
  return(result_project)
}

result_all_project <- lapply(projects, plot_mut_fraction_age, sig_uni = sig_uni, purity = purity)
result_all_project <- do.call(rbind, result_all_project)
result_all_project$q.value <- p.adjust(result_all_project$p.value, method = "BH")
result_all_project$Sig <- ifelse(result_all_project$q.value < 0.05, TRUE, FALSE)
result_all_project$direction <- ifelse(result_all_project$estimate > 0, "up", "down")
write.csv(result_all_project, "Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_sig_multivariate_fraction_all.csv", row.names = FALSE)

write.csv(result_all_project, "Source_Data/Supplementary_Fig_7b.csv", row.names = FALSE)

### Plotting
result_all_project <- read.csv("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_sig_multivariate_fraction_all.csv")
sig_df <- result_all_project[result_all_project$Sig == TRUE,]

dot_break <- c(min(-log10(sig_df$q.value)), median(-log10(sig_df$q.value)), max(-log10(sig_df$q.value)))

# Supplementary Fig 7b
pdf("Analysis_results/Mutations/4.1_Mut_burden_plot_6_classes/Summary_sig_multivariate_fraction_all.pdf", width = 8, height = 4, useDingbats=FALSE) 
p <- ggplot(data = sig_df, aes(x = project, y = class)) +
  geom_point(aes(fill = direction, size = -log10(q.value)), colour="black", pch=21) +
  scale_y_discrete(limits = rev(levels(sig_df$cancer_type))) +
  scale_fill_manual(values=c('#1D91C0','#a50f15')) +
  scale_size_continuous(range = c(3, 6), breaks = dot_break, labels = format(round(dot_break, 1), nsmall = 1)) +
  xlab("Cancer type") +
  ylab("Substitution class") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        #axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9")) +
  guides(size = guide_legend("-log10(adj. p-val)", order = 1))
print(p)
dev.off()

