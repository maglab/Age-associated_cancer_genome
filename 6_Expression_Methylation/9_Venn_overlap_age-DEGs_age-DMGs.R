### Overlap between age-DEGs and age-DMGs
### Fig. 6b and Supplementary Fig. 9

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(VennDiagram)

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
  
  venn.diagram(
    x = list(met_down, met_up, exp_down, exp_up),
    category.names = c("Methy_Down", "Methy_Up", "Expr_Down", "Expr_Up"),
    filename = paste0("Analysis_results/Methylation/3_Overlap_age-DMGs_age-DEGs/Overlap_age-DMGs_age-DEGs_",project,".tiff", collapse = ""),
    output=TRUE,
    
    main = project,
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.cex = 0.45,
    
    # Output features
    imagetype="tiff",
    height = 1200, 
    width = 1200, 
    resolution = 600,
    compression = "none",
    
    # Circles
    lwd = 2,
    lty = 'solid',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = 0.4,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.30,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  
  print(plot)
}

lapply(projects, Venn_plot)



