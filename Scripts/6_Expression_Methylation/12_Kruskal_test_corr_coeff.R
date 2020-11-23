### test distribution of correlation coefficient by Kruskal followed by Dunn's test

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(ggpubr)
library(rstatix)

### projects with > 150 age-DEGs and > 150 age-DMGs
projects <- c("LGG", "BRCA", "UCEC", "ESCA", "KIRP", "OV", "LIHC", "LAML", "SKCM", "PRAD")

### function to test for each cancer
test_kruskal <- function(project){
  
  print(paste0("Working on: ", project))
  ### read data
  df <- read.csv(paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", 
                        project, "_corr_coeff_methy_expr.csv", collapse = ""))
  
  # kruskal test
  res.kruskal <- kruskal_test(Pearson_cor ~ condition, data = df)
  pwc <- dunn_test(Pearson_cor ~ condition, data = df, p.adjust.method = "bonferroni") # pairwise comparisons using dunn test
  pwc$p.adj.format <- formatC(pwc$p.adj, format = "e", digits = 2)  # format p-values to scientific notation
  
  write.csv(as.data.frame(pwc), paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_pwc_dunn_test.csv", collapse = ""), row.names = FALSE)
  
  ### Fig. 6d and Supplementary Fig. 14a
  pdf(paste0("Analysis_results/Methylation/4_Density_plot_corr_coeff_methy_exp/", project, "_kruskal_test.pdf", collapse = ""), width = 6, height = 4.5) 
  p <- ggviolin(df, x = "condition", y = "Pearson_cor", fill = "condition",
                palette = c("#238b45", "#225ea8", "#bd0026", "#feb24c"),
                add = "boxplot", add.params = list(fill = "white")) +
    scale_x_discrete(limits = c("others", "age_DEGs", "age_DMGs", "age_DMGs_DEGs")) +
    ylim(-1.2, 1.2) +
    ylab("Pearson correlation coefficient") +
    ggtitle(project) +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=15,face="bold"),
          legend.position = "none") +
    stat_pvalue_manual(pwc[pwc$group1 == "age_DMGs_DEGs" | pwc$group2 == "age_DMGs_DEGs",], label = "p.adj.format",
                       hide.ns = FALSE, y.position = c(max(df$Pearson_cor) + 0.15, max(df$Pearson_cor), max(df$Pearson_cor) + 0.3)) +
    labs(
      subtitle = get_test_label(res.kruskal)
    )
  
  print(p)
  dev.off()
    
  return(project)
}

lapply(projects, test_kruskal)

