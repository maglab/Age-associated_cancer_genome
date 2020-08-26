### Mutational burden with age
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(ggplot2)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))


# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### function to read mut burden for each cancer type and combine together
get_mut_burden <- function(project){
  mut_burden_df <- read.csv(paste0("Analysis_results/Mutations/Mut_burden/", project, "_mutational_burdens.csv", collapse = ""))
  mut_burden_df$Tumor_Sample_Barcode <- unlist(lapply(as.character(mut_burden_df$Tumor_Sample_Barcode), clean_id))
  mut_burden_df$log_mut <- log10(mut_burden_df$total) # log transform
  mut_burden_df <- mut_burden_df[,c("Tumor_Sample_Barcode", "log_mut")]
  return(mut_burden_df)
}

PANCAN_mut <- lapply(projects, get_mut_burden)
PANCAN_mut <- do.call(rbind, PANCAN_mut)
PANCAN_mut <- merge(PANCAN_mut, clinical, by.x = "Tumor_Sample_Barcode", by.y = "patient")

# linear regression age and mutational burden
lm_fit <- lm(log_mut ~ age + cancer_type + gender + race, data=PANCAN_mut) 
p_value <- formatC(as.numeric(summary(lm_fit)$coefficients[,4][2]), format = "e", digits = 2)
p_value_1 <- as.numeric(summary(lm_fit)$coefficients[,4][2])
r_squared <- round(summary(lm_fit)$adj.r.squared,2)
r_squared_1 <- summary(lm_fit)$r.squared
coeff <- summary(lm_fit)$coefficients[,1][2]

### Supplementary Fig. 5a
pdf("Analysis_results/Mutations/4_Mut_burden_plot/PANCAN_multivariate_mut_burden.pdf", width = 8, height = 6)
my_label <- paste0("adj. R-squared = ", r_squared, "\np = ", p_value)
p <- ggplot(data = PANCAN_mut, aes(x = age, y = log_mut)) + 
  geom_point(aes(color = cancer_type), size = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  ggtitle("PANCAN") +
  xlab("Age at diagnosis") +
  ylab("log10(total mutations)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face="bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  annotate("label", x=-Inf, y = Inf, size = 5,
           label = my_label, hjust=0, vjust=1)
print(p)
dev.off()

result_df <- tidy(lm_fit)

write.csv(result_df, "Analysis_results/Mutations/4_Mut_burden_plot/PANCAN_multivariate_mut_burden.csv", row.names = FALSE)
