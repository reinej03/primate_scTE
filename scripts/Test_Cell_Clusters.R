#Load Files 
barcode1 <- read.delim("GSE127774_ACC_B_barcodes.tsv", header = F, sep = "\t")
genes1 <- data.table::fread("GSE127774_ACC_B_genes.tsv", header = F, data.table = F)
barcode2 <- read.delim("GSE127774_ACC_C_barcodes.tsv", header = F, sep = "\t")
genes2 <- data.table::fread("GSE127774_ACC_C_genes.tsv", header = F, data.table = F)
barcode3 <- read.delim("GSE127774_ACC_H_barcodes.tsv", header = F, sep = "\t")
genes3 <- read.table("GSE127774_ACC_H_genes.tsv", header = F, sep = "\t")
barcode4 <- read.delim("GSE127774_ACC_M_barcodes.tsv", header = F, sep = "\t")
genes4 <- read.table("GSE127774_ACC_M_genes.tsv", header = F, sep = "\t")

library(Matrix)
matrix1 <- readMM("GSE127774_ACC_B_matrix.mtx")
matrix2 <- readMM("GSE127774_ACC_C_matrix.mtx")
matrix3 <- readMM("GSE127774_ACC_H_matrix.mtx")
matrix4 <- readMM("GSE127774_ACC_M_matrix.mtx")


#Add barcodes and genenames to the expression matrixes 
rownames(matrix1) <- genes1$V1
#rownames(matrix2) <- genes2$V1
#rownames(matrix3) <- genes3$V1
#rownames(matrix4) <- genes4$V1
colnames(matrix1) <- barcode1$V1
colnames(matrix2) <- barcode2$V1
colnames(matrix3) <- barcode3$V1
colnames(matrix4) <- barcode4$V1

#create a seurat object for each species 
library(Seurat)
seurat1 <- CreateSeuratObject(counts = matrix1)
seurat2 <- CreateSeuratObject(counts = matrix2)
seurat3 <- CreateSeuratObject(counts = matrix3)
seurat4 <- CreateSeuratObject(counts = matrix4)

#Add the species identifiers for each species 
seurat1@meta.data$orig.ident <- "B"
seurat2@meta.data$orig.ident <- "C"
seurat3@meta.data$orig.ident <- "H"
seurat4@meta.data$orig.ident <- "M"


Idents(seurat1) <- "B"
Idents(seurat2) <- "C"
Idents(seurat3) <- "H"
Idents(seurat4) <- "M"

#normalize each spechies seperately with a scale factor of 10000 (mendioned in paper)
seurat1<- NormalizeData(seurat1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat2<- NormalizeData(seurat2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat3<- NormalizeData(seurat3, normalization.method = "LogNormalize", scale.factor = 10000)
seurat4<- NormalizeData(seurat4, normalization.method = "LogNormalize", scale.factor = 10000)

#find vaiable feature for each spechies seperately with a number of 2000 features (mendioned in paper)
seurat1 <- FindVariableFeatures(seurat1, nfeatures = 2000)
seurat2 <- FindVariableFeatures(seurat2, nfeatures = 2000)
seurat3 <- FindVariableFeatures(seurat3, nfeatures = 2000)
seurat4 <- FindVariableFeatures(seurat4, nfeatures = 2000)

#scate data and runPCA for it, because this is needed to use reduction="rpca" in FindIntegrationAnchors (nor mentioned in paper!!)
seurat1 <- ScaleData(seurat1, verbose = F)
seurat2 <- ScaleData(seurat2, verbose = F)
seurat3 <- ScaleData(seurat3, verbose = F)
seurat4 <- ScaleData(seurat4, verbose = F)

seurat1 <- RunPCA(seurat1, npcs = 30, verbose = F)
seurat2 <- RunPCA(seurat2, npcs = 30, verbose = F)
seurat3 <- RunPCA(seurat3, npcs = 30, verbose = F)
seurat4 <- RunPCA(seurat4, npcs = 30, verbose = F)

#start of the functions listed in their seurat object!

#reduction="rpca" is a different method than in their seurat object!! usually default should be cca 
seurat_anchors <- FindIntegrationAnchors(object.list = list(seurat1, seurat2), dims = 1:30, reduction="rpca") 

#from here on out function calls are exactly like in their seurat object! 
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30) 

seurat_integrated_scaled <- ScaleData(seurat_integrated, verbose = F)

seurat_integrated_PCA <- RunPCA(seurat_integrated_scaled, npcs = 30, verbose = F)

tsne <- RunTSNE(seurat_integrated_PCA, reduction="pca", dims=1:30)

DimPlot(tsne, reduction = "tsne")

#----------------------------------------------------------------------------------------------------------

#check out one of their seurat objects (for AC)

library(Seurat)
premade_seurat <- readRDS("GSE127774_ACC_seurat.rds")
premade_seurat <- UpdateSeuratObject(object = premade_seurat)

DimPlot(premade_seurat, reduction = "tsne")







