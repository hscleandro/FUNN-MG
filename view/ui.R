#Author: Leandro Correa ~@hscleandro

#library(magrittr)
library(shiny)
library(shinyjs)
library(plotly)
library(DT)

library(mongolite)

source(paste0(getwd(),"/support/packages.R"))

CONFIG_PATH <- paste0(getwd(),"/setup.conf")
config_file <- read.table(file = CONFIG_PATH, sep = "=", header = FALSE)
host <- gsub(" ","\\1",as.character(config_file$V2[1]))
port <- gsub(" ","\\1",as.character(config_file$V2[2]))

conect_mongo <- paste0("mongodb://", host, ":", port)

samples = mongo(collection = "samples", db = "manager", url = conect_mongo,
               verbose = FALSE, options = ssl_options())

#source(file.path("server", "helpers.R"))
fluidPage(
  # ==========================================================================
  useShinyjs(),
  
  # add custom JS and CSS
  singleton(
    tags$head(includeScript(file.path('www', 'message-handler.js')),
              includeScript(file.path('www', 'helper-script.js')),
              includeCSS(file.path('www', 'style.css'))
    )
  ),
  
  # enclose the header in its own section to style it nicer
  div(id = "headerSection",
      titlePanel("Results of executions of FUNN-MG tool"),
      
      # author info
      span(
        # author info
        span(
          span("Create by "),
          a("Leandro Corrêa", href = "http://lattes.cnpq.br/6174607285687399"),
          HTML("&bull;"),
          span("Code"),
          a("on GitHub", href = "https://github.com/hscleandro/metagenomics-database"),
          br(),
          
          span("26 de Abril, 2017")
        )
        #HTML("&bull;"),
        #span("Grupo de Genômica ambiental ITV/DS")
        
      )
  ),
  
  # show a loading message initially
  div(
    id = "loadingContent",
    h2("Loading...")
  ),	
  
  # all content goes here, and is hidden initially until the page fully loads
  # ==========================================================================
  hidden(div(id = "allContent",
             
             # sidebar - filters for the data
             sidebarLayout(
               sidebarPanel(
                 # ============================Projeto =======================================	
                 strong(span("Project:")),
                 
                 selectInput(
                   "project", "",
                   c(samples$distinct("project")),
                   selected = "all"),
                 # ============================Amostra =======================================	
                 strong(span("Sample:")),
                 
                 uiOutput("selectSample"),
                 #selectInput(
                 #   "sample", "",
                 #   c(mongo$distinct("sample")),
                 #   selected = "all"),
                 # ===========================Consulta========================================	
                 strong(span("Type of search:")),
                 #span("Tipo de consulta:"),
                 radioButtons(inputId = "query",
                              label = "",
                              choices = c("Pathways" = "path",
                                          "Ko Families" = "ko"
                                          
                                                                                  ),
                              inline = FALSE),
                 #),
                 # ==========================Resultados=======================================	
                 strong(span("Informations:")),
                 HTML("<br><br>"),
                 shinyjs::useShinyjs(),
                 wellPanel(
                   shiny::tags$head(shiny::tags$style(shiny::HTML(
                     "#text {id: text; font-size: 14px; height: 170px; overflow: auto; margin-left: -14px;}"
                   )),
                   shiny::tags$script(src = 'URL.js')),
                   
                   div(id = "text")
                 ),
                 checkboxInput("proteomic", "Proteomic analysis only", FALSE),
                 checkboxInput("subtable", "Information from the subtable", FALSE),
                 #verbatimTextOutput("value")
                 # ==========================================================================
                 #HTML("<br><br><br><br><br><br><br><br><br><br><br><br><br>"),
                 #strong(span("Indicadores de resultado:")),
                 #span("Tipo de consulta:"),
                 HTML("<br><br><br><br><br>"),
                 strong(span("Taxon selected:")),
 
                 uiOutput("selectedControls"),
                 #br(),
                 #textOutput("selectedText"),
                 #selectInput(
                 #   "taxon", "",
                 #   c(" " = "", mongo$distinct("sample")),
                 #   selected = "all"),
                 # ===========================Consulta========================================	
                 strong(span("Tipo de consulta:")),
                 #span("Tipo de consulta:"),
                 radioButtons(inputId = "taxon",
                              label = "",
                              choices = c("Kingdom" = "kingdom",
                                          "Phylum" = "phylum",
                                          "Class" = "class",
                                          "Order" = "ord",
                                          "Family" = "family",
                                          "Genus" = "genre",
                                          "Species" = "species"),
                              inline = FALSE)
                 #),
                 # ============================Função ======================================
                 
                 # ==========================================================================	
                 
                 #strong(span("Indicadores:")),
                
                 
               ),
               # main panel has two tabs - one to show the data, one to plot it
               mainPanel(wellPanel(
                 
                
                   title = "Tabela", id = "tableTab",
                  
                   #downloadButton("downloadData", "Download table"),
                   h3(strong("Click on a table row for more information")
                   ),
                   br(), 
                   DT::dataTableOutput('dataTable'),
                   br(),
                   h4(strong("Information about the line selected above")
                   ),
                   br(), 
                   DT::dataTableOutput('dataSubTable'),
                   br(), br(),
                   h4(strong("Tax distribution from the line selected above")
                   ),
                   br(), 
                   fluidRow(
                     
                     column(6,
                            
                            plotlyOutput("dataPlot",height = "400px")
                     ),
                     column(6,
                            #plotlyOutput("dataPlot",height = "400px")
                            plotlyOutput("dataPlotSample",height = "400px")
                     )
                   ),  
                   br(), br()
                   #DT::dataTableOutput('taxonTable')
               ))
             )
  ))
)

