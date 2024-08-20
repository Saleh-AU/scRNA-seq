## Clean up environment
rm(list = ls())


## Directory
getwd()
setwd("C:\\Users\\salahmed\\Desktop\\RNA Seq")


## Library and packages
install.packages("Seurat")
install.packages("readr")
library(dplyr)
library(Seurat)
library(patchwork)
library(readr)

## Load data
m1<-Read10X(data.dir = "C:\\Users\\salahmed\\Desktop\\RNA Seq\\Pinkie")
m2<-Read10X(data.dir = "C:\\Users\\salahmed\\Desktop\\RNA Seq\\C57B")
dim(m1)
dim(m2)

## Creat a Seurat Object
m1<-CreateSeuratObject(counts = m1, project = "Pinkie", min.cells = 3, min.features = 200)
m2<-CreateSeuratObject(counts = m2, project = "C57B6", min.cells = 3, min.features = 200)

## Mitochondrial percentage 
m1[["percent.mt"]]<-PercentageFeatureSet(m1,pattern = "^MT-")
m2[["percent.mt"]]<-PercentageFeatureSet(m2,pattern = "^MT-")
head(m1@meta.data,10)
head(m2@meta.data,10)
                               
## Visualization with violin plot
VlnPlot(m1,features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
VlnPlot(m2,features = c("nCount_RNA","nFeature_RNA", "percent.mt"))

plot1 <- FeatureScatter(m1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(m1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3<-FeatureScatter(m1,feature1 = "percent.mt", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3
plot2 + plot3

## Selecting cells and mitochondrial percentage (Subsetting)
m1 <- subset(m1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 0.5)
m2 <- subset(m2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 0.5)

VlnPlot(m1,features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
VlnPlot(m2,features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
head(m1@meta.data, 10)

## Normalizing data
m1<-NormalizeData(m1, normalization.method = "LogNormalize", scale.factor = 10000) 
m2<-NormalizeData(m2, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (genes)
m1<-FindVariableFeatures(m1, selection.method = "vst", nfeatures = 2000) 
m2<-FindVariableFeatures(m2, selection.method = "vst", nfeatures = 2000) 
top10m1 <- head(VariableFeatures(m1), 10)
top10m2 <- head(VariableFeatures(m2), 10)

plot1 <- VariableFeaturePlot(m1)
plot2 <- LabelPoints(plot = plot1, points = top10m1, repel = TRUE)
plot1 + plot2
plot3 <- VariableFeaturePlot(m1)
plot4 <- LabelPoints(plot = plot3, points = top10m2, repel = TRUE)
plot2 + plot4

m1m2.list <- list(m1, m2)


# select features that are repeatedly variable across datasets for integration
m1m2.features <- SelectIntegrationFeatures(object.list = m1m2.list) ############


## Performing integration between seurat objects
m1m2.anchors<-FindIntegrationAnchors(object.list = m1m2.list, anchor.features = m1m2.features)
m1m2<-IntegrateData(anchorset=m1m2.anchors)


####### Perform integration analysis
DefaultAssay(m1m2) <- "integrated"

## Scaling data for Linear reduction with PCA
m1m2<-ScaleData(m1m2, verbose = FALSE)

## Linear dimension reduction with PCA
m1m2<-RunPCA(m1m2,npcs = 30,verbose = FALSE)
VizDimLoadings(m1m2,dims = 1:2,reduction = "pca")
DimPlot(m1m2,reduction = "pca")
DimHeatmap(m1m2,reduction = "pca")
DimHeatmap(m1m2,dims = 1:20,cells = 500, reduction = "pca")
ElbowPlot(m1m2)

## Clustering
m1m2<-FindNeighbors(m1m2,dims = 1:15)
m1m2<-FindClusters(m1m2, reductiion = "pca", resolution = 0.8)
head(Idents(m1m2),5)

## Run UMAP or tSNE (Non-linear reduction)
m1m2<-RunUMAP(m1m2,dims = 1:15)
DimPlot(m1m2, reduction = "umap")
p1 <- DimPlot(m1m2, reduction = "umap" , group.by = "orig.ident")
p2 <- DimPlot(m1m2, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2

pdf("cluster.pdf")
DimPlot(m1m2, reduction = "umap",label = TRUE, repel = TRUE,  split.by = "orig.ident")
dev.off()


## Finding specific Cluster markers 
DefaultAssay(m1m2) <- "RNA"
cl1<-FindMarkers(m1m2,ident.1 = 1, min.pct = 0.25)
cl2<-FindMarkers(m1m2,ident.1 = 2, min.pct = 0.25)
head(cl2,10)
cl20<-FindMarkers(m1m2,ident.1 = 20, min.pct = 0.25)
head(cl20,10)

cluster4.markers <- FindMarkers(m1m2, ident.1 = 4, ident.2 = c(1, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
cluster0.markers <- FindMarkers(m1m2, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

######Find markers for every cluster#########
m1m2.markers <- FindAllMarkers(m1m2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cp<-m1m2.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(m1m2.markers, file = "m1m2_markers.csv")
write.csv(cp, file = "cp.csv")
## Violin plot
VlnPlot(m1m2, features = c("Serpinb2", "Retnla", "Siglech","Gm42418"))
VlnPlot(m1m2, features = c("Serpinb2", "Cd209a"))
FeaturePlot(m1m2,features = c("Serpinb2", "Cd209a"), label = TRUE)
## Feature plot
FeaturePlot(m1m2,features = c("Serpinb2", "Retnla", "Siglech"), label = TRUE)
FeaturePlot(m1m2,features = c("Serpinb2", "Cd209a"), label = TRUE)
saveRDS(m1m2,file = "m1m2.rds")

##cells cluster count with identification of markers
m1m2[["my.clusters"]] <- Idents(object = m1m2)

##DEG-across the clusters ((Repetition))###########
m1m2.list <-SplitObject(m1m2, split.by = "orig.ident")
samples <- names(m1m2.list )
master.markers <- NULL
for (i in seq_along(m1m2.list)) {
  tmp.markers <- FindAllMarkers(m1m2.list[[i]], only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  tmp.markers$orig.ident <- samples[i]
  master.markers <- rbind(master.markers, tmp.markers)
}

master.markers %>%
  group_by(cluster, orig.ident) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(tmp.markers, file = "m1m2_DEG_IPA_tmp.csv" )
write.csv(master.markers, file = "m1m2_IPA_master.markers.csv" )############

## Assign cell type ident to clusters
new.cluster.ids <- c("Eosinophil", "Monocyte", "Neutrophil", "Neutrophil", "ILC", "DC", "gd-T cell", "Treg cell", "ILC", "T cell", "NK cell", "Macrophage", "Neutrophil", "DC", "B cell", "Macrophage", "Mast cell", "pre-B", "DC", "Neutrophil", "pDC")
names(new.cluster.ids) <- levels(m1m2)
m1m2 <- RenameIdents(m1m2, new.cluster.ids)
levels(m1m2)
pdf("cell type label.pdf")
DimPlot(m1m2, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5) 
dev.off()
saveRDS(m1m2,file = "m1m2.rds")

## Dot plot
pdf("dot plot.pdf")
Idents(m1m2) <- factor(Idents(m1m2), levels = c("Eosinophil", "Monocyte", "Neutrophil", "Neutrophil", "ILC", "DC", "gd-T cell", "Treg cell", "ILC", "T cell", "NK cell", "Macrophage", "Neutrophil", "DC", "B cell", "Macrophage", "Mast cell", "pre-B", "DC", "Neutrophil", "pDC"))
markers.to.plot <- c("Serpinb2", "Lyz2", "S100a9", "Ccrl2","Il17a","Cd209a","Trdc",
                     "Ctla4", "Il5", "Ccl5","Gzma","Retnla","Il1r2","Cst3", "Igkc", 
                     "Pf4", "Mcpt4", "Stmn1","Ccl22", "Areg", "Siglech")
DotPlot(m1m2, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()

dev.off()

Idents(m1m2) <- factor(Idents(m1m2), levels = c("Eosinophil", "Monocyte", "Neutrophil", "Neutrophil", "ILC", "DC", "gd-T cell", "Treg cell", "ILC", "T cell", "NK cell", "Macrophage", "Neutrophil", "DC", "B cell", "Macrophage", "Mast cell", "pre-B", "DC", "Neutrophil", "pDC"))
markers.to.plot <- c("Serpinb2", "Lyz2", "S100a9", "Ccrl2","Il17a","Cd209a", "Siglech")
DotPlot(m1m2, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()

## Feature plot (Repetition)
pdf("feature plot.pdf")
FeaturePlot(m1m2, features = c("Serpinb2", "Lyz2", "S100a9", "Ccrl2","Il17a","Cd209a","Cd163l1","Ctla4", "Il5", "Ccl5","Gzma","Retnla","Il1r2","Cst3", "Igkc", "Pf4", "Mcpt4", "Stmn1","Ccl22", "S100a8", "Siglech"))
dev.off()
## Violin plot
pdf("vlnplot.pdf")
VlnPlot(m1m2, features = c("Serpinb2", "Lyz2", "S100a9", "Ccrl2","Il17a" ), group.by = "orig.ident",
        pt.size = 1.5, combine = TRUE)
dev.off()

m1m2$celltype <- Idents(m1m2)
pdf("vlnplot2.pdf")
VlnPlot(m1m2, features = c("Il17a"), split.by = "orig.ident", group.by = "celltype",
                 pt.size = 1.2, combine = FALSE)
dev.off()

## Differential gene expression b/n conditions

m1m2<-ScaleData(m1m2)

Idents(m1m2) <- "orig.ident"
strain.markers <- FindMarkers(m1m2, ident.1 = "C57B6", 
                              only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#Selecting Top 40 Genes 
top20 <- strain.markers %>%
  filter(p_val_adj < 0.05) %>% 
  slice_max(abs(avg_log2FC), n = 20) %>% 
  arrange(avg_log2FC) 
top20genes <- rownames(top20)


#Heatmap Across All Cells
DoHeatmap(m1m2, features = top20genes, draw.lines = FALSE)


#Average-based Heatmap
Avg <- AverageExpression(m1m2, return.seurat = TRUE)
DoHeatmap(Avg, features = top40genes, draw.lines = FALSE, slot = "data")




