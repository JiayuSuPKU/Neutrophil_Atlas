## Regulatory gene network analysis
## This file was used to generate input files for pySCENIC

setwd("path")
load("./Score.RData")
load("./preprocessing_Integrated.RData")

library(Matrix)
library(Seurat)

## Load metadata of cells
## Be careful that different versions of metadata might contain different cell names ("ctrl" or "ctl")
cell.info = merge_ref_wt_eco@meta.data
meta.data = read.table("E:/Neutrophil_B/Data/MetaData/ref_meta.txt")

## Randomly select cells for "Stage 1: Estimate co-expression modules"
## Pick non-neutrophils
set.seed(20190719)
NonNeuCells4GRN = c()
for (i in c("B", "DC", "HSPC", "Mono", "T")){
  NonNeuCells4GRN = c(NonNeuCells4GRN, 
                      sample(rownames(meta.data)[meta.data$cell_type==i], 240))
}
table(meta.data[NonNeuCells4GRN,]$cell_type)
NonNeuCells = rownames(meta.data)[meta.data$cell_type %in% 
                                    c("B","DC","HSPC","Mono","T")]

## Pick neutrophils
set.seed(20190719)
Neu4GRN = c()
for (i in levels(factor(cell.info$cluster_condition))){
  Neu4GRN = c(Neu4GRN, sample(
    colnames(merge_ref_wt_eco@data[,cell.info$cluster_condition==i]), 300))
}
table(cell.info[Neu4GRN,]$cluster_condition)

## Load expression profiles from raw data (4 samples)
matrix_dir = "E:/Neutrophil_B/Data/Data1/WTctl_BM/FilteredMatrix/"
mat = Read10X(matrix_dir)
colnames(mat) = paste("wt.ctl.bm1", colnames(mat), sep = "_")
cell.used = colnames(mat) %in% NonNeuCells4GRN
keep.1 = mat[,cell.used]
cell.used = colnames(mat) %in% NonNeuCells 
all.1 = mat[,cell.used]

matrix_dir = "E:/Neutrophil_B/Data/Data2/WTctl_BM2/"
mat = Read10X(matrix_dir)
colnames(mat) = paste("wt.ctl.bm2", colnames(mat), sep = "_")
cell.used = colnames(mat) %in% NonNeuCells4GRN
keep.2 = mat[,cell.used]
cell.used = colnames(mat) %in% NonNeuCells 
all.2 = mat[,cell.used]

matrix_dir = "E:/Neutrophil_B/Data/Data2/WTctl_PB2/FilteredMatrix/"
mat = Read10X(matrix_dir)
colnames(mat) = paste("wt.ctl.pb2", colnames(mat), sep = "_")
cell.used = colnames(mat) %in% NonNeuCells4GRN
keep.3 = mat[,cell.used]
cell.used = colnames(mat) %in% NonNeuCells 
all.3 = mat[,cell.used]

matrix_dir = "E:/Neutrophil_B/Data/Data2/WTctl_SP2/"
mat = Read10X(matrix_dir)
colnames(mat) = paste("wt.ctl.sp2", colnames(mat), sep = "_")
cell.used = colnames(mat) %in% NonNeuCells4GRN
keep.4 = mat[,cell.used]
cell.used = colnames(mat) %in% NonNeuCells 
all.4 = mat[,cell.used]

## Generate expression matrix for network analysis "Stage 1: Estimate co-expression modules"
NonNeuMat4GRN = cbind(keep.1, keep.2, keep.3, keep.4)
genes = intersect(rownames(NonNeuMat4GRN), rownames(merge_ref_wt_eco@raw.data))
Mat4GRN = cbind(NonNeuMat4GRN[genes,], merge_ref_wt_eco@raw.data[genes, Neu4GRN])
## Remove unexpressed genes
genes.filtered = rowSums(Mat4GRN>0)>10
Mat4GRN = Mat4GRN[genes.filtered,]
Cells4GRN = colnames(Mat4GRN)

## Generate expression matrix for network analysis "Stage 3: Calculate AUCell scores for each regulon in each cell"
cell.partition = merge_ref_wt_eco@meta.data$condition == "Control"
mat.ctl = merge_ref_wt_eco@raw.data[, rownames(merge_ref_wt_eco@meta.data)]
mat.ctl = mat.ctl[genes.filtered, cell.partition]
mat.ecoli = merge_ref_wt_eco@raw.data[, rownames(merge_ref_wt_eco@meta.data)]
mat.ecoli = mat.ecoli[genes.filtered, !cell.partition]
mat.nneu = cbind(all.1[genes,], all.2[genes,], 
                 all.3[genes,], all.4[genes,])
mat.nneu = mat.nneu[genes.filtered,]

## Save expression profiles
library(scrattch.io)
write_dgCMatrix_csv(Mat4GRN, "./Mat4GRN.csv", col1_name = "", chunk_size = 1000)
write_dgCMatrix_csv(mat.ctl, "./merge_ref_wt.csv", col1_name = "", chunk_size = 100)
write_dgCMatrix_csv(mat.ecoli, "./merge_ref_ecoli.csv", col1_name = "", chunk_size = 100)
write_dgCMatrix_csv(mat.nneu, "./merge_ref_nneu.csv", col1_name = "", chunk_size = 100)



