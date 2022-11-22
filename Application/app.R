# Load R packages
library(shiny)
library(shinythemes)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(shiny)
library(hdf5r)

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


# Define UI
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
                    tabPanel("UMAP",
                    sidebarPanel(
                      titlePanel("UMAP Input selection"),
                        
                     selectInput(inputId = "cell_type",
                                 label = "Select Event",
                                 choices = c("mouse", "cortex"),
                                 selected = "mouse",
                                 width = "220px"),
                     selectInput(inputId = "Neighbor",
                                 label = "Select Event",
                                 choices = c(10, 20, 30, 40, 50),
                                 selected = 10,
                                 width = "220px"),
                    ),
                    mainPanel(
                      plotOutput(outputId = "UMAP_plot")
                    )
                    ),
                    tabPanel("Interactive Something",
                    "Intentionally left blank",
                      sidebarPanel(
                              tags$h3("Input:"),
                              textInput("txt1", "Given Name:", ""),
                              textInput("txt2", "Surname:", ""),
                              
                            ), # sidebarPanel
                            mainPanel(
                              h1("Header 1"),
                              
                              h4("Output 1"),
                              verbatimTextOutput("txtout"),
                              
                            ), # mainPanel
                  ),
                  ),
                  tabPanel("Hippoccus", "This panel is intentionally left blank"),
                  tabPanel("About us", 
                           fluidRow(
                             shiny::HTML("<br><br><center> 
                                            <h1>About App</h1> 
                                            <h4>What's behind the data.</h4>
                                            </center>"),
                             style = "height:250px;"),
                           fluidRow(
                             div(align = "center",
                                 tags$span(h4("Something"), 
                                           style = "font-weight:bold"
                                 ))
                           ),
                           fluidRow(
                             column(3),
                             column(6,
                                    tags$ul(
                                      tags$li(h6("Bullet points")), 
                                      tags$li(h6("Bullet POintss"))
                                    )
                             ),
                             column(3)
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
      
      brain1 = RunPCA(brain1, verbose = FALSE)
      brain1 <- FindNeighbors(brain1,  reduction = "pca",dims = 1:input$Neighbor)
      brain1 <- FindClusters(brain1, verbose = FALSE)
      brain1 <- RunUMAP(brain1,  reduction = "pca",dims = 1:input$Neighbor)
    }
  )
  output$UMAP_plot <- renderPlot({
    p1 <- DimPlot(brain1, reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(brain1, label = TRUE, label.size = 3)
    p1
  })
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
