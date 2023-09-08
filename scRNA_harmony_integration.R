#devtools::install_github('satijalab/seurat-data')
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

options(timeout = 600) # default timeout is 60 seconds, if your internet is slow you should increase timeout value
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# see available datasets in SeuratData package
AvailableData()

# install the choosen dataset
InstallData("ifnb")

# load dataset to system
LoadData("ifnb")

# get general information related to the dataset
str(ifnb)


# calculate mitochondrial percentage
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)


# filter the dataset 
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)

# data normalization
ifnb.filtered <- NormalizeData(ifnb.filtered)

# finding variable features
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)

# scaling data 
ifnb.filtered <- ScaleData(ifnb.filtered)

# running PCA analysis
ifnb.filtered <- RunPCA(ifnb.filtered)

# observing dimensions with elbow plot
ElbowPlot(ifnb.filtered)

# running UMAP analysis with dimensions from 1 to 20 and according to PCA reduction
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

# creating dimentional plot of results (there is no integration at this step)
before_integration <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')


# start to run harmony package 
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

# see the reductions from harmony object
ifnb.harmony@reductions

# get the harmony embeddings 
ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")

# see harmony embeddings 
ifnb.harmony.embed[1:10,1:10]

# start clustering using UMAP, harmony embeddings will be used instead of PCA embeddings
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# create plot of the clustering results, reduction should be UMAP this will not be problem because it is calculated based on harmony
after_integration <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

# compare the plots 
before_integration|after_integration
