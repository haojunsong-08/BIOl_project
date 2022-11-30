# Load R packages
#install.packages(c("shiny", "shinythemes", "Seurat", "ggplot2", "patchwork", "dplyr", "hdf5r", "visNetwork", "viridis", "ConsensusClusterPlus", "devtools", "ggpubr" , "pheatmap"))
#library(devtools)
#install_github("navinlabcode/CellTrek")
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

#Importing data
data_dir1 = "Anterior/spatial"
data_dir2 = "Anterior" 

#Some other shit loading
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

## Cell_type Distribution
cortex = readRDS('Rds_data/cortex.rds')
cortex.all<-readRDS("Rds_data/allen_cortex.rds")
cortex.all <- SCTransform(cortex.all, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
anchors <- FindTransferAnchors(reference = cortex.all, query = cortex,
                               reference.assay = 'SCT', query.assay = 'SCT',normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = cortex.all$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
## cell_cololization
#setwd(dir="/Users/mosun/Desktop/BIOL 8803/")
cortex_st = readRDS("Rds_data/cortex.rds")

cortex_sc = readRDS("Rds_data/allen_cortex.rds")

cortex_st <- RenameCells(cortex_st, new.names=make.names(Cells(cortex_st)))
cortex_sc <- RenameCells(cortex_sc, new.names=make.names(Cells(cortex_sc)))

#SpatialDimPlot(cortex_st)

brain_traint <- CellTrek::traint(st_data=cortex_st, sc_data=cortex_sc, sc_assay='RNA', cell_names='subclass')
#DimPlot(brain_traint, group.by = "type") 

brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=cortex_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=10, ntree=200, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
#Rshinny visuatlization1
brain_celltrek$subclass <- factor(brain_celltrek$subclass, levels=sort(unique(brain_celltrek$subclass)))

#CellTrek::celltrek_vis(brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, subclass:id_new),
                       #brain_celltrek@images$slice1@image, brain_celltrek@images$slice1@scale.factors$lowres)

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
# Define UI
mst_cons_am <- brain_sgraph_KL_mst_cons
mst_cons_node <- data.frame(id=rownames(mst_cons_am), label=rownames(mst_cons_am))
directed=F  
meta_data= NULL
if (!directed) mst_cons_am[upper.tri(mst_cons_am, diag = T)] <- NA
  
mst_cons_am <- data.frame(id=rownames(mst_cons_am), mst_cons_am, check.names=F)
mst_cons_edge <- reshape2::melt(mst_cons_am) %>% na.omit() %>% magrittr::set_colnames(c('from', 'to', 'value'))

ui <- fluidPage(theme = shinytheme("spacelab"),
                tags$style('.container-fluid {
                             background-color: #f6f3ee;
              }'),
                navbarPage(
                  theme = "spacelab",
                  "ScRNA-Spatial-Mapping",
                  tabPanel("Backgorund", 
                           includeHTML("background.html"),
                           ),
                  navbarMenu("Cortex",
                    tabPanel("Spatial UMAP",
                    sidebarPanel(
                      titlePanel("UMAP Input Selection"),
                        
                     selectInput(inputId = "cell_type",
                                 label = "Tissue",
                                 choices = c("brain", "cortex"),
                                 selected = "mouse",
                                 width = "220px"),
                     selectInput(inputId = "Neighbor",
                                 label = "Dimension",
                                 choices = c(5, 10, 15, 20, 30),
                                 selected = 10,
                                 width = "220px"),
                     textInput(inputId = "gene",
                                 label = "gene_id",
                                 value = "Ttr",
                                 width = "220px"),
                    ),
                    mainPanel(
                      plotOutput(outputId = "UMAP_plot"),
                      plotOutput(outputId = "UMAP_plot2")
                    )
                    ),
                    tabPanel("Single Cell UMAP",
                             sidebarPanel(
                               titlePanel("UMAP Input selection"),
                               
                               selectInput(inputId = "cell_type2",
                                           label = "Select Event",
                                           choices = c("mouse", "cortex"),
                                           selected = "mouse",
                                           width = "220px"),
                               selectInput(inputId = "Neighbor2",
                                           label = "Select Event",
                                           choices = c(5, 10, 15, 20, 30),
                                           selected = 10,
                                           width = "220px"),
                             ),
                             mainPanel(
                               plotOutput(outputId = "UMAP_plot_sc")
                             )
                    ),
                    tabPanel("Cell Type Distribution",
                             sidebarPanel(
                               titlePanel("CellType Feature Selection"),
                               
                               selectInput(inputId = "gene_id_2",
                                           label = "Cell Type",
                                           choices = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", "L6b", "Oligo"),
                                           selected = "Astro",
                                           width = "220px"),
                             ),
                             mainPanel(
                               plotOutput(outputId = "SaptialPlot")
              
                             )
                    ),
                    tabPanel("Cell Colocalization",
                      sidebarPanel(
                        sliderInput('edge_val', 'Edge Value Cutoff',
                                    min=round(range(mst_cons_edge$value)[1], 2), max=round(range(mst_cons_edge$value)[2], 2), value=round(range(mst_cons_edge$value)[1], 2), step=diff(round(range(mst_cons_edge$value), 2))/100),
                        selectInput(inputId='node_col', label='Color', choices=c('None', colnames(brain_cell_class)), selected='None'),
                        selectInput(inputId='node_size', label='Size', choices=c('None', colnames(brain_cell_class)), selected='None'),
                        checkboxInput(inputId='smooth', label='Smooth', value=FALSE),
                        checkboxInput(inputId='physics', label='Physics', value=FALSE),
                        numericInput('mass', 'Mass', value=.5, step=.01),
                        sliderInput('fontsize', 'FontSize', value=15, min=5, max=25),
                        tags$hr(),
                        actionButton('StopID', 'Stop')
                      ),
                      mainPanel(
                        visNetworkOutput("network", height = '800px'))
                  ),
                  tabPanel("Gene Coexpression",
                           sidebarPanel(
                             titlePanel("Cell Type Input selection"),
                             

                             textInput(inputId = "cell",
                                       label = "cell_type",
                                       value = "L5 IT",
                                       width = "220px"),
                             textInput(inputId = "cluster",
                                       label = "cluster_id",
                                       value = "L5 IT VISp ",
                                       width = "220px"),
                           ),
                           mainPanel(
                             plotOutput(outputId = "heatmap")
                           )
                  ),
                  ),
                  tabPanel("About us", 
                           fluidRow(
                             shiny::HTML("<br><br><center> 
                                            <h1>About App</h1> 
                                            <h4>What's behind the Process</h4>
                                            </center>"),
                             style = "height:250px;"),
                           fluidRow(
                             div(align = "center",
                                 tags$img(src = "flowchart.png")
                                 )
                           ),
                           # TEAM BIO
                           fluidRow(
                             column(3),
                             column(6,
                                    shiny::HTML("<br><br><center> <h3>About the team</h3> </center><br>"),
                                    shiny::HTML("<h4>Here is our Team of 4, we are a team of diverse interestand group of students from Georgia Inst
                                    itue of Techonlogy. And here is a little information 
                                                   about the project team!</h4>")
                             ),
                             column(3)
                           ),
                           
                           fluidRow(
                             
                             style = "height:50px;"),
                           
                           fluidRow(
                             column(3),
                             column(2,
                                    div(class="panel panel-default", 
                                        div(class="panel-body",  width = "600px",
                                            align = "center",
                                            div(
                                              tags$img(src = "Mo_sun.jpeg", 
                                                       width = "100px", height = "100px")
                                            ),
                                            div(
                                              tags$h3("Mo Sun"),
                                              tags$h5( tags$i("PHD Bioinformtaics       at Georgia Institute of Technology"))
                                            )
                                        )
                                    )
                             ),
                          
                             column(2,
                                    div(class="panel panel-default", 
                                        div(class="panel-body",  width = "600px",
                                            align = "center",
                                            div(
                                              tags$img(src = "selfie.jpeg", 
                                                       width = "100px", height = "100px")
                                            ),
                                            div(
                                              tags$h3("Haojun Song"),
                                              tags$h5( tags$i("MS CSE and Bioinformtaics       at Georgia Institute of Technology"))
                                            ),
                                            div(
                                              "I love learning and studying new things. "
                                            )
                                        )
                                    )
                             ),
                             
                             
                             column(2,
                                    div(class="panel panel-default", 
                                        div(class="panel-body",  width = "600px",
                                            align = "center",
                                            div(
                                              tags$img(src = "upaasana.jpg", 
                                                       width = "100px", height = "100px")
                                            ),
                                            div(
                                              tags$h3("Upaasana Krishnan"),
                                              tags$h5( tags$i("MS Bioinformtaics       at Georgia Institute of Technology"))
                                            )
                                        )
                                    )
                             ),
                             
                             
                             column(2,
                                    div(class="panel panel-default", 
                                        div(class="panel-body",  width = "600px",
                                            align = "center",
                                            div(
                                              tags$img(src = "nidhi.jpg", 
                                                       width = "100px", height = "100px")
                                            ),
                                            div(
                                              tags$h3("Nidhi Koundinya"),
                                              tags$h5( tags$i("MS Bioinformtaics       at Georgia Institute of Technology"))
                                            )
                                        )
                                    )
                             ),
                             
                           ),
                           fluidRow(style = "height:150px;"))
                  
                ), # navbarPage
                div(style = "margin-bottom: 30px;"),
                tags$footer(
            column(8, "© This is a shinyApp created by class of 2022 at GATECH"),
            column(2, tags$a(href="hsong343@gatech.edu", tags$b("Contact us!"), 
                              class="externallink", style = "color: white; text-decoration: none")), 
       
            column(1, actionLink("twitter_share", label = "star us", icon = icon("github"),
                              style= "color:white;", onclick = sprintf("window.open('%s')", 
                              "https://github.com/haojunsong-08/BIOl_project"))), 
            style = "
            position:fixed;
            text-align:center;
            left: 0;
            bottom:0;
            width:100%;
            z-index:1000;  
            height:30px; /* Height of the footer */
            color: black;
            padding: 10px;
            font-weight: bold;
            background-color: grey"
              ) 
) # fluidPage


# Define server function  
server <- function(input, output) {
  result_UMAP <- reactive(
    {
      brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)
      brain1 = RunPCA(brain1, assay = "SCT", verbose = FALSE)
      brain1 <- FindNeighbors(brain1,  reduction = "pca",dims = 1:input$Neighbor)
      brain1 <- FindClusters(brain1, verbose = FALSE)
      brain1 <- RunUMAP(brain1,  reduction = "pca",dims = 1:input$Neighbor)
    }
  )
  output$UMAP_plot <- renderPlot({
    p1 <- DimPlot(result_UMAP(), reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(result_UMAP(), label = TRUE)
    p1+p2
  })
  output$UMAP_plot2 <- renderPlot({
    SpatialFeaturePlot(result_UMAP(), features = input$gene)
  })
  
  result_UMAP2 <- reactive(
    {
      cortex.all
    })
  output$UMAP_plot_sc <- renderPlot({
    p1 <- DimPlot(result_UMAP2(), group.by = "subclass",label = TRUE, reduction = 'umap')
    p1
  }
  )
  output$SaptialPlot <- renderPlot({
    SpatialFeaturePlot(cortex, features = input$gene_id_2, pt.size.factor = 1.6, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
  })
  
  output$network <- renderVisNetwork({
    if (!is.null(meta_data)) {
      if (input$node_col=='color') {
        col_colmn <- as.character(input$node_col)
        col_df <- data.frame(id=meta_data$id, color=meta_data[, col_colmn])
        mst_cons_node <- dplyr::left_join(mst_cons_node, col_df, by = "id") %>% data.frame
      } else if (input$node_col!='None') {
        col_colmn <- as.character(input$node_col)
        col_df <- data.frame(id=meta_data$id, col_=meta_data[, col_colmn])
        node_cols <- ggpubr::get_palette('Set1', length(unique(col_df$col_)))
        names(node_cols) <- unique(col_df$col_)
        col_df$color <- node_cols[col_df$col_]
        mst_cons_node <- dplyr::left_join(mst_cons_node, col_df, by = "id") %>% data.frame
      }
      if (input$node_size!='None') {
        size_colmn <- as.character(input$node_size)
        size_df <- data.frame(id=meta_data$id, value=as.numeric(meta_data[, size_colmn])/sum(as.numeric(meta_data[, size_colmn])))
        mst_cons_node <- dplyr::left_join(mst_cons_node, size_df, by = "id") %>% data.frame
      }
    }
    mst_cons_edge <- mst_cons_edge[mst_cons_edge$value > input$edge_val, ]
    visNetwork::visNetwork(mst_cons_node, mst_cons_edge) %>%
      visNetwork::visNodes(mass=input$mass, size=15, font = list(size=input$fontsize)) %>%
      visNetwork::visEdges(smooth=input$smooth, physics=input$physics,
                           arrows=list(to=list(enabled=directed, scaleFactor=.2)),
                           color=list(color='rgba(132, 132, 132, .5)'))
  })
  
  output$heatmap <- renderPlot({
    brain_celltrek_l5 <- subset(brain_celltrek, subset=subclass==input$cell)
    brain_celltrek_l5@assays$RNA@scale.data <- matrix(NA, 1, 1)
    brain_celltrek_l5$cluster <- gsub(input$cluster, '', brain_celltrek_l5$cluster)
    
    brain_celltrek_l5 <- FindVariableFeatures(brain_celltrek_l5)
    vst_df <- brain_celltrek_l5@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
    nz_test <- apply(as.matrix(brain_celltrek_l5[['RNA']]@data), 1, function(x) mean(x!=0)*100)
    hz_gene <- names(nz_test)[nz_test<20]
    mt_gene <- grep('^Mt-', rownames(brain_celltrek_l5), value=T)
    rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_l5), value=T)
    vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
    feature_temp <- vst_df$id[1:2000]
    brain_celltrek_l5_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_l5, assay='RNA', approach='cc', gene_select = feature_temp, sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)
    l = list()
    for (n in 1:length(brain_celltrek_l5_scoexp_res_cc$gs)){
      l[[n]] <- rbind(data.frame(gene=c(brain_celltrek_l5_scoexp_res_cc$gs[[n]]), G= n))
    }
    brain_celltrek_l5_k = do.call(rbind,l)%>% 
      magrittr::set_rownames(.$gene) %>% dplyr::select(-1)
    pheatmap::pheatmap(brain_celltrek_l5_scoexp_res_cc$wcor[rownames(brain_celltrek_l5_k), rownames(brain_celltrek_l5_k)], 
                       clustering_method='ward.D2', annotation_row=brain_celltrek_l5_k, show_rownames=F, show_colnames=F, 
                       treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                       color=viridis(10), main= 'spatial co-expression')
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
