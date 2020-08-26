### Plot multiple logistic regression for association between age and arm-level recurrent SCNA results
### Fig. 2d-e
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)

# read gain result
gain_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_gain.csv")
del_df <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_del.csv")

gain_sig <- gain_df[gain_df$Sig == TRUE,]
dim(gain_sig)
gain_sig$p.value <- NULL
gain_sig$direction <- ifelse(gain_sig$estimate > 0, "increase", "decrease")
gain_sig$cancer_type <- as.factor(as.character(gain_sig$cancer_type))
colnames(gain_sig) <- c("cancer_type", "arm", "Regression_coefficient", "std.error",
                        "conf.low", "conf.high", "df.error", "q.value", "Sig", "direction")

del_sig <- del_df[del_df$Sig == TRUE,]
dim(del_sig)
del_sig$p.value <- NULL
del_sig$direction <- ifelse(del_sig$estimate > 0, "increase", "decrease")
del_sig$cancer_type <- as.factor(as.character(del_sig$cancer_type))
colnames(del_sig) <- c("cancer_type", "arm", "Regression_coefficient", "std.error",
                       "conf.low", "conf.high", "df.error", "q.value", "Sig", "direction")


### dot plot for gain
dot_break <- c(min(-log10(gain_sig$q.value)), median(-log10(gain_sig$q.value)), max(-log10(gain_sig$q.value)))

pdf("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_gain.pdf", width = 8, height = 4) 
p <- ggplot(data = gain_sig, aes(x = arm, y = cancer_type)) +
  geom_point(aes(fill = Regression_coefficient, size = -log10(q.value)), colour="black", pch=21) +
  scale_x_discrete("chromosomal arm", limits = c("1q","3p","3q","5p","5q","6p","7p","7q","8q","10p","12p",
                                                 "12q","14q","15q","16p","16q","19p","19q","20p","20q")) +
  scale_y_discrete(limits = rev(levels(gain_sig$cancer_type))) +
  scale_fill_gradient2(mid = "#bdbdbd", midpoint = 0, low = "#1D91C0", high = "#a50f15", guide = "colourbar") +
  scale_size_continuous(range = c(2, 4), breaks = dot_break, labels = format(round(dot_break, 2), nsmall = 2)) +
  ggtitle("Association between age and arm-level copy-number gain") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9")) +
  guides(size = guide_legend("-log10(adj. p-val)", order = 1))
print(p)
dev.off()


### dot plot for del
dot_break <- c(min(-log10(del_sig$q.value)), median(-log10(del_sig$q.value)), max(-log10(del_sig$q.value)))

pdf("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_arm_new/Multivariate_age_recurrent_arm_del.pdf", width = 8, height = 4) 
p <- ggplot(data = del_sig, aes(x = arm, y = cancer_type)) +
  geom_point(aes(fill = Regression_coefficient, size = -log10(q.value)), colour="black", pch=21) +
  scale_x_discrete("chromosomal arm", limits = c("1p","4p","4q","6p","6q","8p","9p","9q","10p","10q","11p",
                                                 "11q","13q","14q","15q","16p","16q","17p","17q","18p","18q","22q")) +
  scale_y_discrete(limits = rev(levels(del_sig$cancer_type))) +
  scale_fill_gradient2(mid = "#bdbdbd", midpoint = 0, low = "#1D91C0", high = "#a50f15", guide = "colourbar") +
  scale_size_continuous(range = c(2, 4), breaks = dot_break, labels = format(round(dot_break, 2), nsmall = 2)) +
  ggtitle("Association between age and arm-level copy-number loss") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_line(colour = "#d9d9d9"),
        panel.grid.minor = element_line(colour = "#d9d9d9")) +
  guides(size = guide_legend("-log10(adj. p-val)", order = 1))
print(p)
dev.off()

