library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/aishw/OneDrive/Desktop/Semester 3/BIOL 8803/Group Project/scRNA-seq Dataset/GSM4829411/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

##Standard pre-processing workflow

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc<-NormalizeData(object=pbmc) 
pbmc<-FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=2000)
pbmc<-ScaleData(object=pbmc)
pbmc<-RunPCA(object=pbmc, dims=1:15)

ElbowPlot(pbmc)

pbmc<-FindNeighbors(object=pbmc, dims=1:15)
pbmc<-FindClusters(object=pbmc, resolution=0.25)
pbmc<-RunUMAP(pbmc, reduction="pca", dims=1:15)
#DimPlot(pbmc, label=T, label.box=T, shuffle=TRUE, raster=FALSE, repel=TRUE, label.size=5)+NoAxes()

#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet", "RBC")
#names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
jpeg('C:/Users/aishw/OneDrive/Desktop/Semester 3/BIOL 8803/Group Project/scRNA-seq Dataset/Results/scRNA-seq_Results.jpeg')
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(pbmc, file = "C:/Users/aishw/OneDrive/Desktop/Semester 3/BIOL 8803/Group Project/scRNA-seq Dataset/Results/scRNA-seq.rds")



