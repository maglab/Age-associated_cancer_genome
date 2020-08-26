### Pathway alterations ComplexHeatmap Cell Cycle in LGG (Fig. 5c)

setwd("/Users/kasitchatsirisupachai/Desktop/Age_differences_cancer/")

library(ggplot2)
library(ComplexHeatmap)

### read data
clinical <- read.csv("Data/all_clin_XML.csv")

### pathway alterations 
pathway_alterations <- read.csv("Data/Pathway_alterations_alteration_level.csv")

clean_id <- function(id){
  tmp <- unlist(strsplit(id, split = "-"))[1:3]
  return(paste0(tmp, collapse = "-"))
}
pathway_alterations$SAMPLE_BARCODE <- unlist(lapply(as.character(pathway_alterations$SAMPLE_BARCODE), clean_id))


### results from multiple regression analysis
df <- read.csv("Analysis_results/Pathway_alterations/Summary_age_multivariate_pathway_alterations.csv")
df <- df[df$Sig == "TRUE",]


# Cell Cycle in LGG
project <- "LGG"
pathway <- "Cell.Cycle"

tmp_pathway_df <- pathway_alterations[pathway_alterations$SAMPLE_BARCODE %in% clinical[clinical$cancer_type == project,]$patient,]
tmp_pathway_df <- tmp_pathway_df[complete.cases(tmp_pathway_df),]
dim(tmp_pathway_df)

CellCycle_genes <- c("CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CCND1", "CCND2",
                     "CCND3", "CCNE1", "CDK2", "CDK4", "CDK6", "RB1", "E2F1", "E2F3")

tmp <- lapply(CellCycle_genes, grepl, colnames(tmp_pathway_df))

TRUE_indices <- unlist(lapply(tmp, which))
TRUE_indices <- unique(sort(TRUE_indices))
tmp_pathway_df <- tmp_pathway_df[,c(1,TRUE_indices)]  # select only columns with alterations in the component of the pathway of interest

rownames(tmp_pathway_df) <- tmp_pathway_df$SAMPLE_BARCODE
tmp_pathway_df$SAMPLE_BARCODE <- NULL
tmp_pathway_df <- tmp_pathway_df[, apply(tmp_pathway_df, 2, sum) > 0]

# function to add alterations to the values in a table
add_alteration <- function(name, tmp_pathway_df){
  alt <- unlist(strsplit(name, split = "[.]"))[1]
  tmp_pathway_df[name] <- ifelse(tmp_pathway_df[name] == TRUE, alt, "")
}
tmp_pathway_df_1 <- lapply(colnames(tmp_pathway_df), add_alteration, tmp_pathway_df = tmp_pathway_df)
tmp_pathway_df_1 <- as.data.frame(do.call(cbind,tmp_pathway_df_1))

# join columns with a same gene
get_gene <- function(name, tmp_pathway_df_1){
  gene <- unlist(strsplit(name, split = "[.]"))[2]
  return(gene)
}
genes <- unlist(lapply(colnames(tmp_pathway_df_1), get_gene))

dup_genes <- genes[duplicated(genes)]   # TP53

join_function <- function(my_row){
  if(length(my_row[my_row != ""]) == 1){
    return(my_row[my_row != ""])
  } else if(length(my_row[my_row != ""]) > 1){
    tmp_row <- my_row[my_row != ""]
    return(paste(tmp_row, collapse = ";"))
  } else if(length(my_row[my_row != ""]) < 1){
    return("")
  }
}

tmp_pathway_df_2 <- tmp_pathway_df_1[,genes == dup_genes[1]]
res_CDKN2A <- apply(tmp_pathway_df_2, 1, join_function)

tmp_pathway_df_2 <- tmp_pathway_df_1[,genes == dup_genes[2]]
res_CDKN2C <- apply(tmp_pathway_df_2, 1, join_function)

tmp_pathway_df_2 <- tmp_pathway_df_1[,genes == dup_genes[3]]
res_RB1 <- apply(tmp_pathway_df_2, 1, join_function)


tmp_pathway_df_1$CDKN2A <- res_CDKN2A
tmp_pathway_df_1$CDKN2C <- res_CDKN2C
tmp_pathway_df_1$RB1 <- res_RB1

tmp_pathway_df_1$DEL.CDKN2A <- NULL
tmp_pathway_df_1$DEL.CDKN2C <- NULL
tmp_pathway_df_1$MUT.CDKN2A <- NULL
tmp_pathway_df_1$MUT.CDKN2C <- NULL
tmp_pathway_df_1$EPISIL.CDKN2A <- NULL
tmp_pathway_df_1$DEL.RB1 <- NULL
tmp_pathway_df_1$MUT.RB1 <- NULL


colnames(tmp_pathway_df_1) <- c("CCND1", "CCND2", "CCND3", "CCNE1", "CDK4", "CDK6", "E2F1",
                                "E2F3", "CDKN2B", "CDKN1B", "CRB1", "CDKN2A", "CDKN2C", "RB1")

tmp_pathway_df_1 <- t(tmp_pathway_df_1)

# annotation df
annot_df <- clinical[clinical$cancer_type == project, c("patient", "age")]
annot_df <- annot_df[annot_df$patient %in% colnames(tmp_pathway_df_1),]
annot_df <- annot_df[order(annot_df$age, decreasing = FALSE),]  # sort by age
age <- annot_df$age

# sort mut_type by age
tmp_pathway_df_1 <- tmp_pathway_df_1[,as.character(annot_df$patient)]
tmp_pathway_df_1 <- as.matrix(tmp_pathway_df_1)

# annotation age and mut burden
ha <- columnAnnotation(age = age,
                       col = list(age = colorRamp2(c(min(age), max(age)), c("#edf8b1", "#1d91c0"))),
                       annotation_legend_param = list(age = list(title = "age", 
                                                                 at = c(min(age), median(age), max(age)),
                                                                 labels = c(min(age), median(age), max(age)))))
ha


col <- c(AMP = "#e31a1c", MUT = "#33a02c", DEL = "#1d91c0", FUSION = "#feb24c", EPISIL = "#5e4fa2")

alter_fun = list(
  AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["AMP"], col = NA)),
  DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["DEL"], col = NA)),
  FUSION = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                          gp = gpar(fill = col["FUSION"], col = NA)),
  MUT = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                       gp = gpar(fill = col["MUT"], col = NA)),
  EPISIL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                          gp = gpar(fill = col["EPISIL"], col = NA))
)

heatmap_legend_param = list(title = "Alterations", at = c("AMP", "DEL", "FUSION", "MUT", "EPISIL"), 
                            labels = c("Gain", "Deletion", "Fusion", "Mutations", "Epigenetic silencing"))

height_pdf = (length(genes) * 0.2) + 2
pdf(paste0("Analysis_results/Pathway_alterations/4_Heatmap_pathway/", project, "_pathway_heatmap_", pathway, ".pdf"), width = 12, height = height_pdf)
p <- oncoPrint(tmp_pathway_df_1, alter_fun = alter_fun, col = col,
               remove_empty_columns = FALSE,
               column_order = as.character(annot_df$patient),top_annotation = ha,
               column_title = paste0(project, ": Cell Cycle", collapse = ""),
               heatmap_legend_param = heatmap_legend_param)
print(p)
dev.off()
