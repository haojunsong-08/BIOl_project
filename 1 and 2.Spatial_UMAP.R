rm(list=ls())
setwd(dir="/Users/mosun/Desktop/BIOL 8803/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(shiny)
library(hdf5r)


data_dir1 = "/Users/mosun/Desktop/BIOL 8803/spatial/Anterior/spatial"
data_dir2 = "/Users/mosun/Desktop/BIOL 8803/spatial/Anterior"


Anterior.brain = Seurat::Read10X_Image(data_dir1, image.name = 'tissue_lowres_image.png')

brain1 = Load10X_Spatial(data.dir=data_dir2,
                         filename = 'V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix.h5',
                         assay='Spatial',
                         slice='slice1',
                         image = Anterior.brain
)

brain1@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(brain1@images[["slice1"]]@coordinates[["tissue"]])
brain1@images[["slice1"]]@coordinates[["row"]] <- as.integer(brain1@images[["slice1"]]@coordinates[["row"]])
brain1@images[["slice1"]]@coordinates[["col"]] <- as.integer(brain1@images[["slice1"]]@coordinates[["col"]])
brain1@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(brain1@images[["slice1"]]@coordinates[["imagerow"]])
brain1@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(brain1@images[["slice1"]]@coordinates[["imagecol"]])
brain1@project.name <-"anterior"
Idents(brain1) <-"anterior"
brain1$orig.ident <-"anterior"


#plot1 <- VlnPlot(brain1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(brain1, features = "nCount_Spatial") + theme(legend.position = "right")
#wrap_plots(plot1, plot2)
# individual 
brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)

brain1 = RunPCA(brain1, verbose = FALSE)
brain1 <- FindNeighbors(brain1, dims = 1:30)
brain1 <- FindClusters(brain1, verbose = FALSE)
brain1 <- RunUMAP(brain1, dims = 1:30)

#brain1 <- FindAllMarkers(brain1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#head(x = brain1)
#dim(brain1)
#data_chr <- filter(brain1, brain1$cluster %in% 0)
#R-shinny(picture)
p1 <- DimPlot(brain1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain1, label = TRUE, label.size = 3)
p1 + p2
#R-shinny(if we can)
SpatialFeaturePlot(brain1, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain1)


##cortex
cortex <- subset(brain1, idents = c(0, 1, 4, 8, 9))
Key(object =cortex@images$slice1)
img<- GetTissueCoordinates(cortex)

ident <- Idents(cortex)
imgId <- data.frame(img,ident)
#head(imgId)
#SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
#p_plot <- ggplot(img, aes(x=imagecol, y=600-imagerow,color = ident)) +
  #geom_point(size = 0.6)
#p_img+p_plot
cortex <- subset(cortex, slice1_imagerow > 400 | slice1_imagecol < 150, invert = TRUE)
#SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
cortex <- subset(cortex, slice1_imagerow > 275 & slice1_imagecol > 370, invert = TRUE)
#SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
cortex <- subset(cortex, slice1_imagerow > 250 & slice1_imagecol > 440, invert = TRUE)

#Rshinny(picture)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2





