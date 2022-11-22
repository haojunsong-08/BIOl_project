setwd(dir="/Users/mosun/Desktop/BIOL 8803/")
cortex = readRDS('cortex.rds')
cortex.all<-readRDS("allen_cortex.rds")
allen_reference <- SCTransform(cortex.all, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

DimPlot(cortex.all, group.by = "subclass",label = TRUE)

anchors <- FindTransferAnchors(reference = cortex.all, query = cortex,
                               reference.assay = 'SCT', query.assay = 'SCT',normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = cortex.all$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1.6, ncol = 2, crop = FALSE, alpha = c(0.1, 1))