
#read in matrices 
library(Seurat)
matrix1 <- Read10X("/Users/joreiner/Desktop/Semester 6/Bachelorarbeit/H")
matrix2 <- Read10X("/Users/joreiner/Desktop/Semester 6/Bachelorarbeit/C")
matrix3 <- Read10X("/Users/joreiner/Desktop/Semester 6/Bachelorarbeit/B")
matrix4 <- Read10X("/Users/joreiner/Desktop/Semester 6/Bachelorarbeit/M")

#name barcodes correctly (spcific to the species)
colnames(matrix1) = paste0("H_",colnames(matrix1))
colnames(matrix2) = paste0("C_",colnames(matrix2))
colnames(matrix3) = paste0("B_",colnames(matrix3))
colnames(matrix4) = paste0("M_",colnames(matrix4))

#filter for just the barcodes that are in the premade seurat object
premade_seurat <- readRDS("premade_seurats/GSE127774_ACC_seurat.rds")
premade_seurat <- UpdateSeuratObject(object = premade_seurat)

barcodes_premade_seurat <- data.frame(colnames(premade_seurat))
colnames(barcodes_premade_seurat) <- c("Barcodes")

barcodes_premade_seurat_H <- data.frame(barcodes_premade_seurat[grepl("^H_", barcodes_premade_seurat$Barcodes), ])
colnames(barcodes_premade_seurat_H) <- c("Barcodes")
barcodes_premade_seurat_C <- data.frame(barcodes_premade_seurat[grepl("^C_", barcodes_premade_seurat$Barcodes), ])
colnames(barcodes_premade_seurat_C) <- c("Barcodes")
barcodes_premade_seurat_B <- data.frame(barcodes_premade_seurat[grepl("^B_", barcodes_premade_seurat$Barcodes), ])
colnames(barcodes_premade_seurat_B) <- c("Barcodes")
barcodes_premade_seurat_M <- data.frame(barcodes_premade_seurat[grepl("^M_", barcodes_premade_seurat$Barcodes), ])
colnames(barcodes_premade_seurat_M) <- c("Barcodes")

select_barcodes_H <-  c(barcodes_premade_seurat_H$Barcodes)
barcode_indices_H <- match(select_barcodes_H, colnames(matrix1))
select_barcodes_C <-  c(barcodes_premade_seurat_C$Barcodes)
barcode_indices_C <- match(select_barcodes_C, colnames(matrix2))
select_barcodes_B <-  c(barcodes_premade_seurat_B$Barcodes)
barcode_indices_B <- match(select_barcodes_B, colnames(matrix3))
select_barcodes_M <-  c(barcodes_premade_seurat_M$Barcodes)
barcode_indices_M <- match(select_barcodes_M, colnames(matrix4))

filtered_matrix_H <- matrix1[,barcode_indices_H]
filtered_matrix_C <- matrix2[,barcode_indices_C]
filtered_matrix_B <- matrix3[,barcode_indices_B]
filtered_matrix_M <- matrix4[,barcode_indices_M]

#create a seurat object for each species 
library(Seurat)
seurat1 <- CreateSeuratObject(counts = filtered_matrix_H)
seurat2 <- CreateSeuratObject(counts = filtered_matrix_C)
seurat3 <- CreateSeuratObject(counts = filtered_matrix_B)
seurat4 <- CreateSeuratObject(counts = filtered_matrix_M)

#normalize each spechies seperately with a scale factor of 10000 
seurat1<- NormalizeData(seurat1, normalization.method = "LogNormalize", scale.factor = 10000)
seurat2<- NormalizeData(seurat2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat3<- NormalizeData(seurat3, normalization.method = "LogNormalize", scale.factor = 10000)
seurat4<- NormalizeData(seurat4, normalization.method = "LogNormalize", scale.factor = 10000)

#find vaiable feature for each spechies seperately with a number of 2000 features 
seurat1 <- FindVariableFeatures(seurat1,selection.method = "vst",  nfeatures = 2000)
seurat2 <- FindVariableFeatures(seurat2,selection.method = "vst" ,nfeatures = 2000)
seurat3 <- FindVariableFeatures(seurat3,selection.method = "vst", nfeatures = 2000)
seurat4 <- FindVariableFeatures(seurat4, selection.method = "vst",nfeatures = 2000)

#find anchors for the integration process  
seurat_anchors2 <- FindIntegrationAnchors(object.list = list(seurat1, seurat2, seurat3, seurat4), dims = 1:30) #reduction="rpca"

#integrate the seurat objects of the diffrent species according to the anchors that were found  
seurat_integrated <- IntegrateData(anchorset = seurat_anchors2, dims = 1:30) 

#scale the integrated data
seurat_integrated <- ScaleData(seurat_integrated, verbose = F)

#perform dimension reduction (PCA) on the integrated data  
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = F)

#run tsne on the integrated integrated data 
seurat_integrated <- RunTSNE(seurat_integrated, reduction="pca", dims=1:30)

#plot what was calculated through running tsne 
DimPlot(seurat_integrated, reduction = "tsne")

#look at marker genes mentioned in paper: 

#create the Feature Plot and save it 
plot <- FeaturePlot(object = seurat_integrated, features = c("GJA1"), pt.size = 1) +
  scale_color_gradientn(colours = c("lightgrey", "blue"), na.value = "lightgrey")

#filter out the rows with NA for the marker gene in qustion 
#plot$data <- na.omit(plot$data)

# sort the data so the dots with a high expression for the gene in question are on top 
plot$data <- plot$data[order(plot$data[[4]], na.last = FALSE),]

#disolay the plot
plot








