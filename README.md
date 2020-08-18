# Steps involved in Subsetting clusture of interest and reclusturing

## Step1 : Loading and plotting the  analysed  data
```
library(Seurat)
CombinedData <- readRDS("CombinedData_labeled.rds")
png(filename = "CombinedData.png")
DimPlot(CombinedData, reduction = "umap")
dev.off()
```

## Step2 : Subsetting the clustures of interest

```
CombinedDataSubset<-subset(CombinedData, ident = c("CD14 Mono", "CD16 Mono", "pDCs", "cDCs"))
```

## Check the Subset plot

```
png(filename = "CombinedDataSubset.png")
DimPlot(CombinedDataSubset, reduction = "umap")
dev.off()
saveRDS(CombinedDataSubset,file="CombinedDataSubset.Rds")
```

## Step3 :Re-clusturing the subset

```
CombinedDataSubset <- ScaleData(CombinedDataSubset, verbose = FALSE)
CombinedDataSubset <- RunPCA(CombinedDataSubset, npcs = 30, verbose = FALSE)
CombinedDataSubset <- RunUMAP(CombinedDataSubset, reduction = "pca", dims = 1:20)
CombinedDataSubset <- RunTSNE(CombinedDataSubset, reduction = "pca", dims = 1:20)
CombinedDataSubset <- FindNeighbors(CombinedDataSubset, reduction = "pca", dims = 1:20)
CombinedDataSubset <- FindClusters(CombinedDataSubset, resolution = 1.0)#resolution was set to 1.0 because this is a large dataset

```
## Step4 : Plot the new clustures based on the subset

```
png(filename = "CombinedData_Reclustured.png")
DimPlot(CombinedDataSubset, reduction = "umap")
dev.off()
```
