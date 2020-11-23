### Overlap between age-DEGs and age-DMGs
### Fig. 6b and Supplementary Fig. 13

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(VennDiagram)
library(grDevices)

### projects with > 150 age-DEGs and > 150 age-DMGs
projects <- c("LGG", "BRCA", "UCEC", "ESCA", "KIRP", "OV", "LIHC", "LAML", "SKCM", "PRAD")

### function to plot venn diagrams
Venn_plot <- function(project){
  
  # read files
  methylation <- read.csv(paste0("Analysis_results/Methylation/1_Methylation_with_age/", project, "_methylation_with_age.csv", collapse = ""))
  met_up <- as.character(methylation[methylation$Sig == TRUE & methylation$estimate > 0,]$gene)
  met_down <- as.character(methylation[methylation$Sig == TRUE & methylation$estimate < 0,]$gene)
  
  expression <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", project, "_gene_expression_with_age.csv", collapse = ""))
  exp_up <- as.character(expression[expression$Sig == TRUE & expression$estimate > 0,]$gene)
  exp_down <- as.character(expression[expression$Sig == TRUE & expression$estimate < 0,]$gene)
  
  myCol <- c("#74c476", "#fed976", "#1d91c0", "#e31a1c")
  
  # write source data
  my_list <- list(met_down, met_up, exp_down, exp_up)
  names(my_list) <- c("Methy_Down", "Methy_Up", "Expr_Down", "Expr_Up")
  my_list <- melt(my_list)
  colnames(my_list) <- c("Gene", "Category")
  
  if(project %in% c("LGG", "BRCA")){
    write.csv(my_list, paste0("Source_Data/Fig_6b_", project, ".csv"), row.names = FALSE)
  } else {
    write.csv(my_list, paste0("Source_Data/Supplementary_Fig_13_", project, ".csv"), row.names = FALSE)
  }
  
  head(my_list)
  myplot <- venn.diagram(
    x = list(met_down, met_up, exp_down, exp_up),
    category.names = c("Methy_Down", "Methy_Up", "Expr_Down", "Expr_Up"),
    #filename = paste0("Analysis_results/Methylation/3_Overlap_age-DMGs_age-DEGs/Overlap_age-DMGs_age-DEGs_",project,".tiff", collapse = ""),
    filename = NULL, 
    output=TRUE,
    
    main = project,
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 1.4,
    
    # Output features
    imagetype="tiff",
    height = 300, 
    width = 300, 
    resolution = 600,
    compression = "none",
    
    # Circles
    lwd = 2,
    lty = 'solid',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = 1.3,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.7,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  
  pdf(paste0(paste0("Analysis_results/Methylation/3_Overlap_age-DMGs_age-DEGs/Overlap_age-DMGs_age-DEGs_",project,".pdf", collapse = "")), width = 5, height = 5, useDingbats=FALSE)
  grid.draw(myplot)
  dev.off()
}

lapply(projects, Venn_plot)



