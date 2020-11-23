### Oncoprint of age biases in mutations (multiple logistic regrssion)
### Fig. 4g and Supplementary Fig. 9
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(ComplexHeatmap)
library(ggplot2)
library(circlize)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")

### sig mutations
sig_mut <- read.csv("Analysis_results/Mutations/Summary_age_SNVs_multivariate_new.csv")
sig_mut <- sig_mut[sig_mut$Sig == TRUE,]
sig_mut$direction <- ifelse(sig_mut$estimate > 0, "increase", "decrease")

projects <- sort(unique(as.character(sig_mut$cancer_type)))

# function to convert id
convert_id <- function(id){
  tmp <- unlist(strsplit(id, split = "[.]"))
  return(paste0(tmp, collapse = "-"))
}
# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to create oncoprint for each project

mut_oncoprint <- function(project, clinical){
  
  print(paste("Start working on: ", project))
  
  df <- sig_mut[sig_mut$cancer_type == project,]
  genes <- as.character(df$gene)
  
  ### mut burden
  mut_burden <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", 
                                project, "_mutational_burdens.csv", collapse = ""))
  mut_burden$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden$Tumor_Sample_Barcode), clean_id))
  
  ### clinical
  clin <- clinical[clinical$cancer_type == project,]
  
  ### SNV type table
  mut_type <- read.csv(paste0("Analysis_results/Mutations/3.2_Age_SNVs_mut_types/", 
                              project, "_mutation_types_new.csv", collapse = ""))
  mut_type[is.na(mut_type)] <- ""
  rownames(mut_type) <- as.character(mut_type$Hugo_Symbol)
  mut_type$Hugo_Symbol <- NULL
  colnames(mut_type) <- unlist(lapply(as.character(colnames(mut_type)), convert_id))
  mut_type <- mut_type[rownames(mut_type) %in% genes,]
  
  # annotation df
  annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
  annot_df <- annot_df[annot_df$patient %in% colnames(mut_type),]
  annot_df <- merge(annot_df, mut_burden, by.x = "patient", by.y = "Tumor_Sample_Barcode")
  annot_df$log_mut <- log10(annot_df$total)
  annot_df <- annot_df[order(annot_df$age, decreasing = FALSE),]  # sort by age
  age <- annot_df$age
  total_mut <- as.numeric(annot_df$log_mut)
  
  # sort mut_type by age
  mut_type <- mut_type[,as.character(annot_df$patient)]
  mut_type <- as.matrix(mut_type)
  
  # write source data
  if(project %in% c("LGG", "GBM")){
    write.csv(mut_type, paste0("Source_Data/Fig_4g_", project, ".csv"))
  } else {
    write.csv(mut_type, paste0("Source_Data/Supplementary_Fig_9_", project, ".csv"))
  }
  
  # annotation age and mut burden
  ha <- columnAnnotation(age = age, 
                         mutational_burden = anno_points(total_mut, ylim = c(0, ceiling(max(total_mut))), axis = TRUE),
                         col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                         annotation_legend_param = list(age = list(title = "age", 
                                                                   at = c(min(age), median(age), max(age)),
                                                                   labels = c(min(age), median(age), max(age)),
                                                                   title_gp = gpar(fontsize=13, fontface="bold"),
                                                                   labels_gp = gpar(fontsize=13))))
  ha
  
  # annotation increase or decrease freq with age
  annot_direction <- df[,c("gene", "direction")]
  
  row_ha <- rowAnnotation(direction = as.character(annot_direction$direction),
                          col = list(direction = c("increase" = "#fdae61", "decrease" = "#7bccc4")),
                          annotation_legend_param = list(title_gp = gpar(fontsize=13, fontface="bold"),
                                                         labels_gp = gpar(fontsize=13)))
  
  col <- c(Frame_Shift_Del = "#6a3d9a", Frame_Shift_Ins = "#1f78b4", In_Frame_Del = "#fb9a99",
           In_Frame_Ins = "#80b1d3", Missense_Mutation = "#33a02c", Nonsense_Mutation = "#e31a1c",
           Nonstop_Mutation = "#feb24c", Splice_Site = "#ff7f00", Translation_Start_Site = "#cc4c02",
           Multiple = "#000000")
  
  alter_fun <- list(
    Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                     gp = gpar(fill = col["Frame_Shift_Del"], col = NA)),
    Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                     gp = gpar(fill = col["Frame_Shift_Ins"], col = NA)),
    In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["In_Frame_Del"], col = NA)),
    In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["In_Frame_Ins"], col = NA)),
    Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["Missense_Mutation"], col = NA)),
    Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["Nonsense_Mutation"], col = NA)),
    Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                      gp = gpar(fill = col["Nonstop_Mutation"], col = NA)),
    Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["Splice_Site"], col = NA)),
    Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                            gp = gpar(fill = col["Translation_Start_Site"], col = NA)),
    Multiple = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["Multiple"], col = NA)))
  
  height_pdf = (length(genes) * 0.5) + 2
  pdf(paste0("Analysis_results/Mutations/3.3_Age_SNVs_heatmap/", project, "_age_SNVs_heatmap_new.pdf"), width = 12, height = height_pdf, useDingbats=FALSE) 
  p <- oncoPrint(mut_type, alter_fun = alter_fun, col = col,
                 remove_empty_columns = FALSE,
                 right_annotation = row_ha,
                 column_order = as.character(annot_df$patient),
                 top_annotation = ha,
                 column_title = project,
                 column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                 row_title_gp = gpar(fontsize = 14),
                 row_names_gp = gpar(fontsize = 14),
                 pct_gp = gpar(fontsize = 14),
                 heatmap_legend_param = list(title = "Alterations",
                                             title_gp = gpar(fontsize = 13, fontface = "bold"), 
                                             labels_gp = gpar(fontsize = 13)))
  print(p)
  dev.off()

}

lapply(projects, mut_oncoprint, clinical = clinical)


