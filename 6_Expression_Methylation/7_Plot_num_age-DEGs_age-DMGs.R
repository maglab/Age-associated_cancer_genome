### plot number of age-DEGs and age-DMGs

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(reshape)

### read age-DEGs and age-DMGs
expression <- read.csv("Analysis_results/Gene_Expression/Summary_gene_expression_with_age.csv")
methylation <- read.csv("Analysis_results/Methylation/Summary_DNA_methylation_with_age.csv")

expression <- expression[order(expression$all, decreasing = TRUE),]
expression$all <- NULL
expression <- expression[expression$project != "READ",] # remove READ because samples in methylation < 100
projects <- as.character(expression$project)
expression <- melt(expression)
colnames(expression) <- c("project", "direction", "number_of_genes")
expression[expression$number_of_genes == 0,]$number_of_genes <- NA

methylation$all <- NULL
methylation <- melt(methylation)
colnames(methylation) <- c("project", "direction", "number_of_genes")
methylation[methylation$number_of_genes == 0,]$number_of_genes <- NA


### Fig. 6a
### dot plot for expression
pdf("Analysis_results/Gene_Expression/Number_of_age-DEGs.pdf", width = 8, height = 3) 
p <- ggplot(data = expression, aes(x = project, y = direction)) +
  geom_point(aes(fill = direction, size = number_of_genes), colour="black", pch=21) +
  scale_size_continuous(range = c(1, 6), breaks = c(50, 500, 1000, 2000, 2500)) +
  scale_x_discrete(limits = projects) +
  scale_fill_manual(values = c("#1D91C0", "#a50f15"), labels=c("Down-regulated", "Up-regulated")) +
  ggtitle("Number of age-DEGs") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9"))
print(p)
dev.off()

### dot plot for methylation
pdf("Analysis_results/methylation/Number_of_age-DMGs.pdf", width = 8, height = 3) 
p <- ggplot(data = methylation, aes(x = project, y = direction)) +
  geom_point(aes(fill = direction, size = number_of_genes), colour="black", pch=21) +
  scale_size_continuous(range = c(1, 6), breaks = c(100, 500, 1000, 2000, 3000)) +
  scale_x_discrete(limits = projects) +
  scale_fill_manual(values = c("#1D91C0", "#a50f15"), labels=c("Down-regulated", "Up-regulated")) +
  ggtitle("Number of age-DMGs") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 9),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9"))
print(p)
dev.off()

