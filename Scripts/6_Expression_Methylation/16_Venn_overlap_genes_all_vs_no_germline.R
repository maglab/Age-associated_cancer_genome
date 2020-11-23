### Overlap between age-DEGs identified from all samples and from samples excluding samples with germline variants

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(VennDiagram)
library(ggplot2)
library(grDevices)

### Venn diagram of overlapping genes
Venn_plot <- function(project, sig_level_all, sig_level_no_germline){
  
  # read result from all samples
  All <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  All_Sig <- as.character(All[All$q.value < sig_level_all,]$gene)
  
  # read result from samples without germline mutations
  no_germline <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age_exclude_germline.csv", collapse = ""))
  no_germline_Sig <- as.character(no_germline[no_germline$q.value < sig_level_no_germline,]$gene)
  
  myCol <- c("#2c7fb8", "#edf8b1")
  
  if(length(All_Sig) >= length(no_germline_Sig)){
    my_rotation = 0
  } else {
    my_rotation = 180
  }
  
  if(sig_level_no_germline == 0.05){
    my_label = "No germline"
  } else if (sig_level_no_germline == 0.1) {
    my_label = "No germline (p<0.1)"
  }
  
  # write source data
  my_list <- list(All_Sig, no_germline_Sig)
  names(my_list) <- c("All_Sig", "no_germline_Sig")
  my_list <- melt(my_list)
  colnames(my_list) <- c("Gene", "Category")
  
  write.csv(my_list, paste0("Source_Data/Supplementary_Fig_12b_expression_", project, ".csv", collapse = ""), row.names = FALSE)
  
  myplot <- venn.diagram(
    x = list(All_Sig, no_germline_Sig),
    category.names = c("All samples", my_label),
    #filename = paste0("Analysis_results/Gene_Expression/Overlap_age-DEGs_all_and_exclude_germline_",project, "_", sig_level_all, "_", sig_level_no_germline, ".tiff", collapse = ""),
    filename = NULL, 
    output=TRUE,
    
    main = paste0(project, ": Expression", collapse = ""),
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 1.4,
    
    # Output features
    imagetype="tiff",
    height = 300, 
    width = 300, 
    resolution = 600,
    
    # Circles
    lwd = 2,
    lty = 'solid',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = 1.4,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.pos = c(0, 0),
    
    scaled = FALSE,
    rotation.degree = my_rotation)
  
  pdf(paste0("Analysis_results/Gene_Expression/Overlap_age-DEGs_all_and_exclude_germline_",project, "_", sig_level_all, "_", sig_level_no_germline, ".pdf", collapse = ""),
      width = 4, height = 4, useDingbats = FALSE)
  grid.draw(myplot)
  dev.off()
}

lapply(c("BRCA", "OV", "UCEC"), Venn_plot, sig_level_all  = 0.05, sig_level_no_germline = 0.05)


### Overlap between age-DMGs identified from all samples and from samples excluding samples with germline variants

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(VennDiagram)
library(ggplot2)

### Venn diagram of overlapping genes
Venn_plot <- function(project, sig_level_all, sig_level_no_germline){
  
  # read result from all samples
  All <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  All_Sig <- as.character(All[All$q.value < sig_level_all,]$gene)
  
  # read result from samples without germline mutations
  no_germline <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age_exclude_germline.csv", collapse = ""))
  no_germline_Sig <- as.character(no_germline[no_germline$q.value < sig_level_no_germline,]$gene)
  
  myCol <- c("#2c7fb8", "#edf8b1")
  
  if(length(All_Sig) >= length(no_germline_Sig)){
    my_rotation = 0
  } else {
    my_rotation = 180
  }
  
  if(sig_level_no_germline == 0.05){
    my_label = "No germline"
  } else if (sig_level_no_germline == 0.1) {
    my_label = "No germline (p<0.1)"
  }
  
  # write source data
  my_list <- list(All_Sig, no_germline_Sig)
  names(my_list) <- c("All_Sig", "no_germline_Sig")
  my_list <- melt(my_list)
  colnames(my_list) <- c("Gene", "Category")
  
  write.csv(my_list, paste0("Source_Data/Supplementary_Fig_12b_methylation_", project, ".csv", collapse = ""), row.names = FALSE)
  
  myplot <- venn.diagram(
    x = list(All_Sig, no_germline_Sig),
    category.names = c("All samples", my_label),
    #filename = paste0("Analysis_results/Methylation/Overlap_age-DMGs_all_and_exclude_germline_",project, "_", sig_level_all, "_", sig_level_no_germline, ".tiff", collapse = ""),
    filename = NULL,
    output=TRUE,
    
    main = paste0(project, ": Methylation", collapse = ""),
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 1.4,
    
    # Output features
    imagetype="tiff",
    height = 300, 
    width = 300, 
    resolution = 600,
    
    # Circles
    lwd = 2,
    lty = 'solid',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = 1.4,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.pos = c(0, 0),
    
    scaled = FALSE,
    rotation.degree = my_rotation)
  pdf(paste0("Analysis_results/Methylation/Overlap_age-DMGs_all_and_exclude_germline_",project, "_", sig_level_all, "_", sig_level_no_germline, ".pdf", collapse = ""),
      width = 4, height = 4, useDingbats=FALSE)
  grid.draw(myplot)
  dev.off()
}

lapply(c("BRCA", "OV", "UCEC"), Venn_plot, sig_level_all  = 0.05, sig_level_no_germline = 0.05)


