library(shiny)
library(shinythemes)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(visNetwork)
library(viridis)
library(ConsensusClusterPlus)
library(CellTrek)
setwd(dir = "/Users/haojunsong/Downloads/BIOLProejct//BIOL_project/Application/Rds_data")
cortex = readRDS('cortex.rds')
cortex.all<-readRDS("allen_cortex.rds")
cortex.all <- SCTransform(cortex.all, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

#DimPlot(cortex.all, group.by = "subclass",label = TRUE)

anchors <- FindTransferAnchors(reference = cortex.all, query = cortex,
                               reference.assay = 'SCT', query.assay = 'SCT',normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = cortex.all$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = "Astro", pt.size.factor = 1.6, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

