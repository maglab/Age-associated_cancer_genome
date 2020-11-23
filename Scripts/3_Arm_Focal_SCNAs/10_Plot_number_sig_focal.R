### plot number of sig focal regions
### Fig. 3a
setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

library(ggplot2)
library(ggrepel)
library(reshape2)

# read significant results
gain <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_gain_sig.csv")

loss <- read.csv("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/Multivariate_age_recurrent_focal_loss_sig.csv")

### plot number of regions
# gain
gain$direction <- ifelse(gain$estimate > 0, "increase", "decrease")
gain_df <- gain[,c("cancer_type", "direction")]
gain_df <- table(gain_df)
gain_df <- melt(gain_df)

write.csv(gain_df, "Source_Data/Fig_3a_gain.csv", row.names = FALSE)

pdf("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/No_of_age_recurrent_focal_gain_sig.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = cancer_type, y = value, fill = direction), data = gain_df) +
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title="Gain", y = "Number of regions") +
  scale_fill_manual(values=c('#1D91C0','#a50f15')) +
  scale_x_discrete(sort(unique(as.character(gain_df$cancer_type)))) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()


# loss
loss$direction <- ifelse(loss$estimate > 0, "increase", "decrease")
loss_df <- loss[,c("cancer_type", "direction")]
loss_df <- table(loss_df)
loss_df <- melt(loss_df)

write.csv(loss_df, "Source_Data/Fig_3a_loss.csv", row.names = FALSE)

pdf("Analysis_results/CNAs/3_Age_biases_recurrent_SCNAs/Age_recurrent_focal_new/No_of_age_recurrent_focal_loss_sig.pdf", width = 6, height = 4.5, useDingbats=FALSE) 
p <- ggplot(aes(x = cancer_type, y = value, fill = direction), data = loss_df) +
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title="Loss", y = "Number of regions") +
  scale_fill_manual(values=c('#1D91C0','#a50f15')) +
  scale_x_discrete(sort(unique(as.character(gain_df$cancer_type)))) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()
