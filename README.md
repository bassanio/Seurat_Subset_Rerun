# Steps involved in Subsetting clusture of interest and reclusturing

[Method1](#Method1)

[Method2](#Method2)

# Subset Clustures
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


# Method1 


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

# Method2

Continuation after step2 from Method 1
## Step3 : Change teh defualt assay
```
DefaultAssay(CombinedDataSubset) <- "RNA"
```

## Step4: Split the Data for reanalysis based on Visits

For downstream analysis each sample/condition need minimum of 200 cells. else its better to remove those samples or merge with other samples
```
Subset_Cells.list <- SplitObject(CombinedDataSubset, split.by = "Visits")
```
## Step5: Rerun SCT
```
for (i in 1:length(Subset_Cells.list)) {
  Subset_Cells.list[[i]] <- PercentageFeatureSet(Subset_Cells.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  Subset_Cells.list[[i]] <- SCTransform(Subset_Cells.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
  #Subset_Cells.list[[i]] <- NormalizeData(Subset_Cells.list[[i]], verbose = FALSE)
  #Subset_Cells.list[[i]] <- FindVariableFeatures(Subset_Cells.list[[i]], selection.method = "vst", nfeatures = 2000,  verbose = FALSE)
}
```

## Step6:finding anchors and standard processing 
```
pag_combined.anchors <- FindIntegrationAnchors(object.list = Subset_Cells.list, dims = 1:20)
CombinedDataSubset.combined <- IntegrateData(anchorset = pag_combined.anchors, dims = 1:20)
CombinedDataSubset.combined <- ScaleData(CombinedDataSubset.combined, verbose = FALSE)
CombinedDataSubset.combined <- RunPCA(CombinedDataSubset.combined, npcs = 6, verbose = FALSE)
CombinedDataSubset.combined <- RunUMAP(CombinedDataSubset.combined, reduction = "pca", dims = 1:6)
CombinedDataSubset.combined <- RunTSNE(CombinedDataSubset.combined, reduction = "pca", dims = 1:6,check_duplicates = FALSE)
CombinedDataSubset.combined <- FindNeighbors(CombinedDataSubset.combined, reduction = "pca", dims = 1:6)
CombinedDataSubset.combined <- FindClusters(CombinedDataSubset.combined)
```
## Step7: plot the UMAP
```
png(filename = "CombinedData_Reclustured_BYVisit.png")
DimPlot(CombinedDataSubset.combined, reduction = "umap")
dev.off()
```

