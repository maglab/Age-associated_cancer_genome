### Check number of up- and down-regulated genes with age

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer")

# clinical file
clinical <- read.csv("Data/all_clin_XML.csv")
projects <- unique(as.character(clinical$cancer_type))

# remove projects that have samples < 100
projects <- projects[!(projects %in% c("ACC", "CHOL", "DLBC", "KICH", "MESO", "READ", "THYM", "UCS", "UVM"))]

check_age_genes <- function(project){
  df <- read.csv(paste0("Analysis_results/Gene_Expression/1_Gene_expression_with_age/", 
                        project, "_gene_expression_with_age.csv", collapse = ""))
  df <- df[df$Sig == TRUE,]
  down <- nrow(df[df$estimate < 0,])
  up <- nrow(df[df$estimate > 0,])
  
  out <- cbind(project, down, up)
  return(out)
}

results <- lapply(projects, check_age_genes)
results <- do.call(rbind, results)
results <- as.data.frame(results)
results$project <- as.character(results$project)
results$down <- as.numeric(as.character(results$down))
results$up <- as.numeric(as.character(results$up))
results$all <- results$down + results$up

write.csv(results, "Analysis_results/Gene_Expression/Summary_gene_expression_with_age.csv", row.names = FALSE)
