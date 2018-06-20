#Author: Leandro Correa ~@hscleandro

library(shiny)
#library(magrittr)
#library(plyr)
library(dplyr)
#library(tidyr)
library(ggplot2)
library(shinyjs)
library(plotly)
library(DT)
library(jsonlite)
library(mongolite)
library(devtools)

source(paste0(getwd(),"/support/packages.R"))

CONFIG_PATH <- paste0(getwd(),"/setup.conf")
config_file <- read.table(file = CONFIG_PATH, sep = "=", header = FALSE)
host <- gsub(" ","\\1",as.character(config_file$V2[1]))
port <- gsub(" ","\\1",as.character(config_file$V2[2]))

conect_mongo <- paste0("mongodb://", host, ":", port)

samples = mongo(collection = "samples", db = "manager", url = conect_mongo,
              verbose = FALSE, options = ssl_options())

#local = mongo(collection = "sequences", db = "local", url = conect_mongo,
#             verbose = FALSE, options = ssl_options())
#kos_seq <- local$find(query = '{"$and":[{"kegg_ko":{"$ne":"nan", "$exists": "true"}, "phylum" : "Verrucomicrobia"}]}',
#                      field = '{"kegg_ko":""}'
#                  )
#kos_seq <- kos_seq$kegg_ko
#kos_seq <- toJSON(kos_seq)
#
#funn2 = mongo(collection = "funn", db = "local", url = conect_mongo,
#              verbose = FALSE, options = ssl_options())

#result <- funn2$find(query = paste0('{"kegg_ko":{"$in":',kos_seq,'}}'))

shinyServer(function(input, output, session) {
  
# ============= BUILDING THE INPUTS ======================================================================================
  
  printInShinyBox <- function(phrase){
    foo <- function(message_teste) {
      message(message_teste)
      Sys.sleep(0.5)
    }
    
    withCallingHandlers({
      shinyjs::html("text", "")
      foo(phrase)
    },
    message = function(m) {
      shinyjs::html(id = "text", html = m$message, add = TRUE)
    })
  }
  
  ### This will create the dynamic dropdown list ###
  
  output$selectSample <- renderUI({
    
    project <- as.character(input$project) 
    #project <- "CANGA"
    sample <- samples$distinct("sample_name", query = paste0('{"project" : "',project,'"}'))
    if(length(sample) == 0){
      selectInput("sample", "", c("",""))  
    }
    else{
      selectInput("sample", "", sample)  
    
    }
    
  })
  
  
  output$selectedControls <- renderUI({
    if(input$subtable == TRUE){
      data <- subDat()
      s <- input$dataSubTable_rows_selected
     
    }
    else{
      data <- cDatTable()
      s <- input$dataTable_rows_selected
     
    }
    
    if(is.null(s)){
      selectInput("taxons", "", c("",""))
      
    }
    else{
      
      data_taxon <- taxonDat()
      selectInput("taxons", "", data_taxon$taxon)
  
    }
    
    
  })
  

# ============= MANIPULATE THE DATA ======================================================================================

	path_Dat <- reactive({
	  sample <- as.character(input$sample)
	  project <- as.character(input$project) 
	  
	  #sample <- "MG_34_EMMA"
	  #project <- "CANGA"
	  
	  collection_name <- paste0(project,"-",sample)
	  
	  funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
	                  verbose = FALSE, options = ssl_options())
	  #sample <- "MG34_MGRAST_N"
	  #project <- "CANGA"
    data <- funn$find(query = paste0('{"pathways": {"$exists": "true", "$ne":"NA"}}'),
	                     field = '{"pathways" : "","class" : "","subclass" : "","Function" : "","genes_relation_sample" : "",
                                 "load_roume" : "","load_rahman" : "","choke_point" : "","p_value":""}')
    
	  ordenation <- order(data$load_roume, decreasing = TRUE)
	  data <- data[ordenation, c("_id","pathways","class","subclass","Function","p_value","genes_relation_sample","load_roume","load_rahman","choke_point")]
	  data$genes_relation_sample <- as.numeric(as.character(data$genes_relation_sample))
	  data$choke_point <- round(as.numeric(as.character(data$choke_point)), digits = 2)
	  data$`_id` <- seq(1:nrow(data))
	  rownames(data) <- NULL
	  
	  data
	})
  
  ko_Dat <- reactive({

    sample <- as.character(input$sample)
    project <- as.character(input$project) 
    #sample <- "MG_34_EMMA"
    #project <- "CANGA"
    
    collection_name <- paste0(project,"-",sample)
    
    funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    kaas = mongo(collection = collection_name, db = "kaas", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    prot = mongo(collection = collection_name, db = "proteomics", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    
    data <- funn$find(query = paste0('{"kegg_ko": {"$exists": "true", "$ne":"NA"}}'),
                       field = '{"kegg_ko" : "","ko_functions" : "","degree" : "","load_roume" : "","load_rahman" : "","choke_point" : ""}')
    
    ordenation <- order(data$load_roume, decreasing = TRUE)
    data <- data[ordenation, c("_id","kegg_ko","ko_functions","degree","load_roume","load_rahman","choke_point")]
    data$degree <- as.numeric(as.character(data$degree))
    data$choke_point <- round(as.numeric(as.character(data$choke_point)), digits = 2)
    data$`_id` <- seq(1:nrow(data))
    rownames(data) <- NULL
    
    samples$find(query = paste0('{"kegg_ko": {"$exists": "true", "$ne":"NA"}}'))
    
    proteomic <- samples$find(query = paste0('{"$and":[{"sample_name": "', sample,'",
    	                                                  "project": "', project,'"}]}'),
                              field = '{"proteomic_tool":""}')
    
    if(input$proteomic){
      
      kaas_ko <- kaas$find(query = paste0('{"kegg_ko": {"$ne":"NA", "$ne":"nan"}}'))
      
      id_seq <- kaas_ko$id_seq
      id_seq <- toJSON(id_seq)
      
      kos_prot <- prot$find(query = paste0('{"id_seq":{"$in":',id_seq,'}}'))
      
      index <- which(kaas_ko$id_seq %in% kos_prot$id_seq)
      
      kos_prot_kegg_ko <- kaas_ko$kegg_ko[index]
      
      data <- data[data$kegg_ko %in% unique(kos_prot_kegg_ko),]
    }
    
    
    
    data
  })
  
  cDatTable <- reactive({
		
    if (input$query == "ko"){
      data <- ko_Dat()
    }
    if (input$query == "path"){
      data <- path_Dat()
    }
		
		#message <- paste0("<b>Tipo de dado: </b>", class(data),"</br>")
		#message <- paste0(message,"</br>","<b>Maior valor captado</b></br> R$",format(max(data$Captacao), nsmall = 2))
		#printInShinyBox(message)
		
	  data
	})
  
  ko_subDat <- reactive({
    data <- cDatTable()
    
    if(is.null(input$dataTable_rows_selected)) {s <- 1}
    else {s = input$dataTable_rows_selected}
    
    s = input$dataTable_rows_selected
    ko <- data$kegg_ko[s]
    
    #ko <- "K01915"
    
    sample <- as.character(input$sample)
    project <- as.character(input$project)
    
    collection_name <- paste0(project,"-",sample)
    
    funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    
    
    id <- funn$find(query = paste0('{"kegg_ko":"',ko,'"}'),
                    field = '{"paths":""}')
    
    paths <- id$paths[[1]]
    data <- path_Dat()
    index_paths <- which(data$pathways %in% paths)
    data <- data[index_paths,]
    
    
    
    data
  })
  
  path_subDat <- reactive({
    data <- cDatTable()
    
    if(is.null(input$dataTable_rows_selected)) {s <- 1}
    else {s = input$dataTable_rows_selected}
    
    path <- data$pathways[s]
    
    #path <- "path:map02020"
    
    sample <- as.character(input$sample)
    project <- as.character(input$project)
    
    collection_name <- paste0(project,"-",sample)
    
    funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
   
    id <- funn$find(query = paste0('{"pathways":"',path,'"}'),
                     field = '{"kos" :""}')
    
    kos <- id$kos[[1]]
    data <- ko_Dat()
    index_kos <- which(data$kegg_ko %in% kos)
    data <- data[index_kos,]
    
    data
  })
  
  subDat <- reactive({
    if (input$query == "ko"){
      data <- ko_subDat()
      
    }
    if (input$query == "path"){
      data <- path_subDat()
    }
    
    data
    
  })
  
  taxonDat <- reactive({

    if(input$subtable == TRUE){
      data <- subDat()
      s <- input$dataSubTable_rows_selected
    }
    else{
      data <- cDatTable()
      s = input$dataTable_rows_selected
    }
    #cat("aqui2\n")
    sample <- as.character(input$sample)
    project <- as.character(input$project) 
    
    collection_name <- paste0(project,"-",sample)
    
    kaas = mongo(collection = collection_name, db = "kaas", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    taxon = mongo(collection = collection_name, db = "kaiju", url = conect_mongo,
                  verbose = FALSE, options = ssl_options())
    funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
                 verbose = FALSE, options = ssl_options())
    
    # message <- paste0("<b>features: </b>", input$subtable,"</br>")
    #  printInShinyBox(message)
    
    if("kegg_ko" %in% colnames(data)){
      #cat("aqui3\n")
      ko <- data$kegg_ko[s]
      #ko <- "K02259" 
      kaas_ko <- kaas$find(query = paste0('{"kegg_ko": "', ko,'"}'))
      
      id_seq <- kaas_ko$id_seq
      id_seq <- toJSON(id_seq)
      
      taxon_info <- taxon$find(query = paste0('{"id_seq":{"$in":',id_seq,'}}'),
                               field = paste0('{"id_seq" :"","kegg_ko" : "","kingdom":"","phylum":"","class" : "","order" : "","family" : "","genre" : "","species" : ""}')
      )
      
    }
    else{

      path <- data$pathways[s]
      #path <- "path:map00195"
      
      id <- funn$find(query = paste0('{"pathways":"',path,'"}'),
                       field = '{"kos" :""}')
      ko <- id$kos[[1]]
      ko <- toJSON(ko)
      #project <- "CANGA"
      #sample <- "MG201"
      kaas_ko <- kaas$find(query = paste0('{"kegg_ko":{"$in":',ko,'}}'))
      
      id_seq <- kaas_ko$id_seq
      id_seq <- toJSON(id_seq)
      
      taxon_info <- taxon$find(query = paste0('{"id_seq":{"$in":',id_seq,'}}'),
                               field = paste0('{"id_seq" :"","kegg_ko" : "","kingdom":"","phylum":"","class" : "","order" : "","family" : "","genre" : "","species" : ""}')
      )
      ko <- path
      
    }
    
    #taxon_info
    
    if(input$taxon == "kingdom"){
      data_taxon <- taxon_info %>%
        group_by(kingdom) %>%
        select(kingdom) %>%
        summarise(freqSequence = n())
      
    }
    else if(input$taxon == "phylum"){
      data_taxon <- taxon_info %>%
        group_by(phylum) %>%
        select(phylum) %>%
        summarise(freqSequence = n())
      
    }
    else if(input$taxon == "class"){
      data_taxon <- taxon_info %>%
        group_by(class) %>%
        select(class) %>%
        summarise(freqSequence = n())
    }
    else if(input$taxon == "ord"){
      data_taxon <- taxon_info %>%
        group_by(order) %>%
        select(order) %>%
        summarise(freqSequence = n())
    }
    else if(input$taxon == "family"){
      data_taxon <- taxon_info %>%
        group_by(family) %>%
        select(family) %>%
        summarise(freqSequence = n())
    }
    else if(input$taxon == "genre"){
      data_taxon <- taxon_info %>%
        group_by(genre) %>%
        select(genre) %>%
        summarise(freqSequence = n())
    }
    else if(input$taxon == "species"){
      data_taxon <- taxon_info %>%
        group_by(species) %>%
        select(species) %>%
        summarise(freqSequence = n())
    }
    
    data_taxon <- as.data.frame(data_taxon)
    colnames(data_taxon) <- c("taxon","freqSequence")
    #index <- which(data_taxon$taxon == "NA")
    
    #if(index > 0){ 
    #  data_taxon <- data_taxon[-index,]
    #  cat(index)
    #}
    
    data_taxon
    
  })
  
# ============= TAB TO SHOW DATA IN TABLE ================================================================================

  output$dataTable <- DT::renderDataTable(
    
    cDatTable(), options = list(lengthChange = FALSE, responsive = TRUE,
                 pageLength = 10), rownames= FALSE, selection = 'single'
  )
  
  output$dataSubTable <- DT::renderDataTable(
   
      subDat(), options = list(lengthChange = FALSE,
                                  pageLength = 5, sDom  = '<"top">lrt<"bottom">ip'), rownames= FALSE, selection = 'single'
    
  )
  
  output$taxonTable <- DT::renderDataTable(
    
    taxonDat(), options = list(lengthChange = FALSE,
                             pageLength = 5, sDom  = '<"top">lrt<"bottom">ip'), rownames= FALSE, selection = 'single'
    
  )
  

# ============= TAB TO PLOT DATA =========================================================================================
	
	buildPlot <- reactive({
  	  
	  project <- as.character(input$project) 
	  sample <- as.character(input$sample)
	  
	  collection_name <- paste0(project,"-",sample)
	  
	  funn = mongo(collection = collection_name, db = "funn", url = conect_mongo,
	               verbose = FALSE, options = ssl_options())
	  prot = mongo(collection = collection_name, db = "proteomics", url = conect_mongo,
	               verbose = FALSE, options = ssl_options())
	  kaas = mongo(collection = collection_name, db = "kaas", url = conect_mongo,
	               verbose = FALSE, options = ssl_options())
	  txn = mongo(collection = collection_name, db = "kaiju", url = conect_mongo,
	               verbose = FALSE, options = ssl_options())
	  
	    if(input$subtable == TRUE){
	      data <- subDat()
	      s <- input$dataSubTable_rows_selected
	    }
	    else{
	      data <- cDatTable()
	      s <- input$dataTable_rows_selected
	    }
	    
    	if(is.null(s)){
	      p <- plotly_empty()
	    }
    	else{
    	  #cat("aqui1\n")
    	  data_taxon <- taxonDat()
    	  #cat("aqui2\n")
    	  if("kegg_ko" %in% colnames(data)){
    	    name <- data$kegg_ko[s]
    	    #name <- "K02259"
    	    kaas_ko <- kaas$find(query = paste0('{"kegg_ko" : "', name,'"}'))
    	    count <- nrow(kaas_ko)
    	    kegg_ko_rest <- paste0("http://www.genome.jp/dbget-bin/www_bget?", name)
    	    
    	    address <- paste0('<a href="',kegg_ko_rest,'">', name,'</a>')
    	    message <- paste0("<b>Grupo ortólogo identificado: </b>", address,"</br>")
    	    message <- paste0(message,"<b>Quantidade total na amostra: </b>", count,"</br>")
    	    #cat("count",count,"\n")
    	    
    	    id_seq <- kaas_ko$id_seq
    	    id_seq <- toJSON(id_seq)
    	    
    	    count <- prot$count(query = paste0('{"id_seq":{"$in":',id_seq,'}}'))
    	    
    	    if(count > 0) {proteomic <- TRUE}
    	    else {proteomic <- FALSE}
    	    message <- paste0(message,"<b>Análise Proteômica: </b>", proteomic,"</br>")
    	    #message <- paste0(message,"<b>Acesso ao banco de dados KEGG: </b>", address,"</br>")
    	  
    	  }
    	  else{
    	    name <- data$pathways[s]
    	    #name <- "path:map00920"
    	    pathways <-funn$find(query = paste0('{"pathways":"', name,'"}'),
    	                           field = '{"pathways" : "","genes_relation_sample" : "","genes_relation_noted" : ""}')
    	    
    	    kegg_ko_rest <- paste0("http://www.genome.jp/dbget-bin/www_bget?", name)
    	    address <- paste0('<a href="',kegg_ko_rest,'">', name,'</a>')
    	    message <- paste0("<b>Metabolic pathway identified: </b>", address,"</br>")
    	    message <- paste0(message,"<b>Number of related groups in the sample: </b>", pathways$genes_relation_sample,"</br>")
    	    message <- paste0(message,"<b>Number of groups annotated in KEGG database: </b>", pathways$genes_relation_noted,"</br>")
    	    
    	  }
    	  #project <- as.character(input$project) 
    	  p <- plot_ly(data_taxon, 
    	               labels = data_taxon$taxon, 
    	               values = data_taxon$freqSequence, 
    	               type = "pie",
    	               source = "PiePlot"
    	  ) %>%
    	  layout(title = paste0("\t\t\tDistribution of taxons - ",name,"          "), showlegend = TRUE)
    	  #xaxis = list(title = "", side = "top"), 
    	  #yaxis = list(title = ""), 
    	  #       autosize = T,
    	  #       width =  400, 
    	  #       height = 200 
    	  #height = resize(length(index)), 
    	  #margin = m
    	  # )
    	  
    	  #message <- paste0(message, "</br>","<b>parcial: </b>", rel,"</br>")
    	  #sample <- "MG34_MGRAST_N"
    	  #project <- "CANGA"
    	  ######################################################################################################################
    	  #total <- sequence$count(query = paste0('{"$and":[{"',input$taxon,'": {"$exists": "true", "$ne":"NA"},
        #                                                "project": "', project,'",
        #                                                "id_sample": "', sample,'"}]}'))
    	  type_taxon <- as.character(input$taxon)
    	  #type_taxon <- "phylum"
    	  
    	  total_tx <- txn$find(query = paste0('{"',type_taxon,'": {"$exists": "true", "$ne":"NA"}}'),
    	                        field = paste0('{"',type_taxon,'" :""}'))
    	  
        total <- length(total_tx[,type_taxon])
        total_tx <- length(unique(total_tx[,type_taxon]))
    	  
    	  data_taxon <- taxonDat()
    	  tx <- data_taxon$taxon
    	  index <- which(tx == "NA")
    	  if(length(index)){ tx <- tx[-index] }
    	  total_rel <- length(unique(tx))
    	  tx <- toJSON(tx)
    	  rel <- txn$count(query = paste0('{"', type_taxon,'":{"$in":',tx,'}}'))
    	  
    	  percet_taxon <- (rel *100) / total
    	  percet_taxon <- round(percet_taxon, digits = 2)
    	  message <- paste0(message,"<b>Percentage of taxons involved: </b>", percet_taxon,"% </br>")
    	  message <- paste0(message,"<b>Number of taxons involved: </b>", total_rel,"/",total_tx,"</br>")
    	  #######################################################################################################################
    	  
    	  printInShinyBox(message)
    	}
    	p

	})
  
  buildPlotSample <- reactive({
    
    project <- as.character(input$project) 
    sample <- as.character(input$sample)
    
    collection_name <- paste0(project,"-",sample)
    
    tx = mongo(collection = collection_name, db = "kaiju", url = conect_mongo,
                  verbose = FALSE, options = ssl_options())
    
    if(input$subtable == TRUE){
      s <- input$dataSubTable_rows_selected
    }
    else{
      s <- input$dataTable_rows_selected
    }
    
    if(is.null(s)){
      p <- plotly_empty()
    }
    else{
      taxon <- as.character(input$taxons)
      type_taxon <- as.character(input$taxon)
      
      total <- tx$count(query = paste0('{"',type_taxon,'": {"$ne":"NA"}}'))
      rel <- tx$count(query = paste0('{"',type_taxon,'" : "', taxon,'"}'))
      
      name <- c("outher taxons",taxon)
      amount <- c(total-rel,rel)
      
      qnt_taxon <- cbind(name, amount)
      qnt_taxon <- as.data.frame(qnt_taxon)
      
      p <- plot_ly(qnt_taxon, 
                   labels = name, 
                   values = amount, 
                   type = "pie",
                   source = "PiePlot"
      ) %>%
      layout(title = paste0("\t\t\tDistribution of the selected taxon in the sample."), showlegend = TRUE)
    }

    p

  })
	
 	output$dataPlot <- renderPlotly({

 	    buildPlot()

	 
	})
 	
 	output$dataPlotSample <- renderPlotly({
 	  
 	  buildPlotSample()
 	  
 	})

# ========== LOADING THE APP ==========
	dataValues <- reactiveValues(
		appLoaded = FALSE
	)
	
	# Show form content and hide loading message
	
	hide("loadingContent")
	show("allContent")
})
