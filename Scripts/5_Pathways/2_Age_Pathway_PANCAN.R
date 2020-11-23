### Pan-cancer age-associated pathway alterations

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(broom)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(scales)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# clean id function
clean_id <- function(id){
  tmp <- gsub(pattern = "[.]", replacement = "-", id)
  tmp <- unlist(strsplit(tmp, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}

### pathway alteration file
pathway_df <- read.csv("Data/Pathway_alterations_clean.csv")
rownames(pathway_df) <- pathway_df$SAMPLE_BARCODE
pathway_df$SAMPLE_BARCODE <- NULL

###################################################################################################
### Logistic regression to test whether age associates with an alteration in each pathway

# function to test association between age and pathway alterations
test_age_pathway <- function(pathway, pathway_df, clinical){
  
  print(paste("Working on: ", pathway))
  # get value for pathway of interest
  values_pathway <- pathway_df[pathway]
  values_pathway$patient <- rownames(values_pathway)
  rownames(values_pathway) <- NULL
  colnames(values_pathway) <- c("alteration", "patient")
  values_pathway <- values_pathway[complete.cases(values_pathway),]
  
  # merge with age data
  df <- merge(values_pathway, clinical, by.x = "patient", by.y = "patient")
  
  # logistic regression
  logit_fit <- glm(alteration ~ age + gender + race + cancer_type , data = df, family = "binomial")
  
  summary(logit_fit)
  
  result <- tidy(logit_fit)
  CI <- confint.default(logit_fit, level = 0.95)
  result <- cbind(as.data.frame(result), CI)
  result <- result[,c("term", "estimate", "std.error", "statistic", "2.5 %", "97.5 %", "p.value")]
  colnames(result) <- c("term", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")
  write.csv(result, paste0("Analysis_results/Pathway_alterations/1_PANCAN_pathway_alterations/PANCAN_pathway_alterations_", pathway, "_glm.csv", collapse = ""), row.names = FALSE)
  
  result_df <- as.data.frame(result[result$term == "age",])
  result_df$term <- pathway
  colnames(result_df) <- c("pathway", "estimate", "std.error", "statistic", "conf.low", "conf.high", "p.value")
  
  return(result_df)
}


pathways <- colnames(pathway_df)

tmp <- do.call(rbind,lapply(pathways, test_age_pathway, pathway_df = pathway_df, clinical = clinical))
tmp$p.value <- as.numeric(as.character(tmp$p.value))
tmp$q.value <- p.adjust(tmp$p.value, method = "BH")
tmp$Sig <- ifelse(tmp$q.value < 0.05, TRUE, FALSE)

tmp$odds <- exp(tmp$estimate)
tmp$odds_conf.low <- exp(tmp$conf.low)
tmp$odds_conf.high <- exp(tmp$conf.high)

tmp <- tmp[,c("pathway", "estimate", "std.error", "conf.low", "conf.high", "statistic",
              "odds", "odds_conf.low", "odds_conf.high", "p.value", "q.value", "Sig")]

write.csv(tmp, "Analysis_results/Pathway_alterations/PANCAN_pathway_alterations_multivariate_summary.csv", row.names = FALSE)


### plot
df <- read.csv("Analysis_results/Pathway_alterations/PANCAN_pathway_alterations_multivariate_summary.csv")
df

### add color
my_colour <- c()
for(i in 1:nrow(df)){
  if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
    my_colour <- c(my_colour, "Red")
  } else {
    my_colour <- c(my_colour, "Grey")
  }
}

df$my_colour <- factor(my_colour)
df$pathway <- c("Cell Cycle", "HIPPO", "MYC", "NOTCH", "NRF2", "PI3K", "RTK/RAS", "TP53", "TGF-Beta", "WNT")

# write source data
write.csv(df, "Source_Data/Fig_5a.csv", row.names = FALSE)

### Fig. 5a
pdf("Analysis_results/Pathway_alterations/PANCAN_pathway_alterations_multivariate_summary.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = estimate, y = -log10(q.value), color = my_colour, label = pathway), data = df) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("#bdbdbd", "#a50f15")) +
  xlab("Regression coefficient") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("PANCAN: Age and pathway alteration") +
  xlim(c(-0.013,0.013)) +
  geom_label_repel(size = 5) +
  geom_hline(
    yintercept = c(-log10(0.05),-log10(0.05)),
    col = "#bdbdbd",
    linetype = "dashed",
    size = 0.5) +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15,face="bold"),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()
