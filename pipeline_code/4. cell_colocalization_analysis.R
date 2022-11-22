setwd(dir="/Users/mosun/Desktop/BIOL 8803/")
library (dplyr)
library (ggplot2)
library (harmony)
library (patchwork)
library (Seurat)
library(CellTrek)
library(dplyr)
library(Seurat)
library(viridis)
library(ConsensusClusterPlus)
library(shiny)

cortex_st = readRDS("cortex.rds")

cortex_sc = readRDS("allen_cortex.rds")

cortex_st <- RenameCells(cortex_st, new.names=make.names(Cells(cortex_st)))
cortex_sc <- RenameCells(cortex_sc, new.names=make.names(Cells(cortex_sc)))

#SpatialDimPlot(cortex_st)

brain_traint <- CellTrek::traint(st_data=cortex_st, sc_data=cortex_sc, sc_assay='RNA', cell_names='subclass')
#DimPlot(brain_traint, group.by = "type") 

brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=cortex_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=30, ntree=500, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
#Rshinny visuatlization1
brain_celltrek$subclass <- factor(brain_celltrek$subclass, levels=sort(unique(brain_celltrek$subclass)))

CellTrek::celltrek_vis(brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, subclass:id_new),
                       brain_celltrek@images$slice1@image, brain_celltrek@images$slice1@scale.factors$lowres)

#Cell colocalization analysis
glut_cell <- c('Astro','CR','Endo','Lamp5','Macrophage','Meis2','L2/3 IT', 'L4', 'L5 IT', 'L5 PT', 'NP', 'L6 IT', 'L6 CT',  'L6b',
               'Oligo','Peri','Pvalb','Serpinf1','SMC','Sncg','Sst','Vip','VLMC')
names(glut_cell) <- make.names(glut_cell)
brain_celltrek_glut <- subset(brain_celltrek, subset=subclass %in% glut_cell)
brain_celltrek_glut$subclass <- factor(brain_celltrek_glut$subclass, levels=glut_cell)
brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, col_cell='subclass', use_method='KL', eps=1e-50)

brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- glut_cell[colnames(brain_sgraph_KL_mst_cons)]
## We then extract the metadata (including cell types and their frequencies)
brain_cell_class <- brain_celltrek@meta.data %>% dplyr::select(id=subclass) %>% unique
brain_celltrek_count <- data.frame(freq = table(brain_celltrek$subclass))
brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")