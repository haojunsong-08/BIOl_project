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

brain_celltrek$subclass <- factor(brain_celltrek$subclass, levels=sort(unique(brain_celltrek$subclass)))



## gene coexpression analysis
#'L5 IT'
brain_celltrek_l5 <- subset(brain_celltrek, subset=subclass=='L5 IT')
brain_celltrek_l5@assays$RNA@scale.data <- matrix(NA, 1, 1)
brain_celltrek_l5$cluster <- gsub('L5 IT VISp ', '', brain_celltrek_l5$cluster)
#DimPlot(brain_celltrek_l5, group.by = "cluster", label = TRUE)

brain_celltrek_l5 <- FindVariableFeatures(brain_celltrek_l5)
vst_df <- brain_celltrek_l5@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_l5[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(brain_celltrek_l5), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_l5), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]
brain_celltrek_l5_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_l5, assay='RNA', approach='cc', gene_select = feature_temp, sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)
brain_celltrek_l5_k <- rbind(data.frame(gene=c(brain_celltrek_l5_scoexp_res_cc$gs[[1]]), G='K1'), 
                             data.frame(gene=c(brain_celltrek_l5_scoexp_res_cc$gs[[2]]), G='K2')) %>% 
  magrittr::set_rownames(.$gene) %>% dplyr::select(-1)

#Rshinny visuatlization3
pheatmap::pheatmap(brain_celltrek_l5_scoexp_res_cc$wcor[rownames(brain_celltrek_l5_k), rownames(brain_celltrek_l5_k)], 
                   clustering_method='ward.D2', annotation_row=brain_celltrek_l5_k, show_rownames=F, show_colnames=F, 
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='L5 IT spatial co-expression')

