library(ggplot2)
library(Matrix)
library(Seurat)
library(copykat)
options(future.globals.maxSize = 200 * 1024^3)
library(future)
plan("multicore", workers = 32)
setwd("/home/chenliqun/copykat/Cle")
scRNA<-readRDS("/home/chenliqun/copykat/Cle/ClearCellRenalCellCarcinoma_singleRCelltype.rds")
scRNA@meta.data$cell_type[scRNA@meta.data$cell_type == "Epithelial_cells"] <- "Epithelial cells (malignant)"
target_cells <- subset(scRNA, subset = cell_type %in% c("Epithelial cells (malignant)"))
target_cells[["RNA"]] <- JoinLayers(target_cells[["RNA"]])
exp.rawdata <- LayerData(target_cells, assay = "RNA", layer = "counts")
raw_counts <- as.matrix(exp.rawdata)
dim(raw_counts)
copykat.test <- copykat(
  rawmat = raw_counts,
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = "/home/chenliqun/copykat/Cle/ClearCellRenalCellCarcinoma",
  distance = "euclidean",
  norm.cell.names = "",
  plot.genes = FALSE,
  genome = "hg20",
  n.cores = 32
)
if (!dir.exists("CopyKAT_plots")) {
  dir.create("CopyKAT_plots")
}
tumor_cells <- colnames(copykat.test$CNAmat)[copykat.test$prediction == "aneuploid"]
pred.test <- data.frame(copykat.test$prediction)
table(pred.test$copykat.pred)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]
CNA.test <- data.frame(copykat.test$CNAmat)
scRNA$CopyKAT = copykat.test$prediction$copykat.pred
tumor_malignant_cells <- subset(scRNA, subset = cell_type %in% "Epithelial cells (malignant)" & CopyKAT == "aneuploid")
expr_matrix <- GetAssayData(tumor_malignant_cells, assay = "SCT", slot = "data")
if (!dir.exists("CopyKAT_Data")) {
  dir.create("CopyKAT_Data")
}
write.csv(expr_matrix,file = "/home/chenliqun/copykat/Cle/CopyKAT_Data/ClearCellRenalCellCarcinoma_all_malignant_cells_expression.csv",row.names = TRUE,quote = FALSE)
if(ncol(tumor_malignant_cells) > 0){
  tumor_expr_matrix <- GetAssayData(tumor_malignant_cells, assay = "SCT", slot = "data")
  write.csv(tumor_expr_matrix,file = "/home/chenliqun/copykat/Cle/CopyKAT_Data/ClearCellRenalCellCarcinoma_malignant_aneuploid_cells_expression.csv",row.names = TRUE,quote = FALSE)
  tumor_expr_matrix_t <- as.data.frame(as.matrix(t(tumor_expr_matrix)))
  write.csv(tumor_expr_matrix_t,file = "/home/chenliqun/copykat/Cle/CopyKAT_Data/ClearCellRenalCellCarcinoma_malignant_aneuploid_cells_expression_transposed.csv",row.names = TRUE,quote = FALSE)
  print(paste("Successfully saved", ncol(tumor_malignant_cells), "aneuploid tumor cell expression matrix"))
} else {
  warning("No aneuploid tumor cells detected!")
}