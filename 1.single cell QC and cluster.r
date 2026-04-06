##Conduct single cell quality control, dimensionality reduction and unsupervised clustering.
library(Seurat)
data <- Read10X(data.dir = 'input') #The variable "input" represents the directory path where the matrix files are stored.
data <- CreateSeuratObject(counts = data ,project = 'seurat', min.cells = 3, min.features = 200,names.delim = '_')
data[["percent.mt"]]<-PercentageFeatureSet(object = data, pattern = "^MT-")
data<-subset(x=data,subset=nFeature_RNA>200&nFeature_RNA<4000&percent.mt<50&nCount_RNA>1000&nCount_RNA<20000) 
obj_list<-SplitObject(input, split.by = "patient")
obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = obj_list)
obj.anchors <- FindIntegrationAnchors(object.list =obj_list, anchor.features = features)
obj.combined <- IntegrateData(anchorset = obj.anchors)
##
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindClusters(obj.combined, resolution = 0.3)