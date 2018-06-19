#===================================================================================================
#
# FUNN-MG: FuNctioNal Network Analysis of MetaGenomics Data
#
# Funn-MG v1.11 (2016-10-11) -- "C'est la vie"
# Copyright 2016 Leandro Corrêa
#
# This file is part of FUNN-MG.
#
# FUNN-MG is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FUNN-MG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#===================================================================================================
#load("/home/leandro/R/version/funn-mg-v.13/data/tempfiles/funn-CANGA-MG_34_EMMA.RData")
#output_status <- FALSE
#source(paste0(getwd(),"/supports/bibfunction.R"))
#library(igraph)

## PT[Br] Escrevendo as informações das vias selecionadas em um arquivo .csv
ord <- order(table_information$load_roume, decreasing = TRUE)
table_information <- table_information[ord,]

## PT[Br] Impriindo um piechart com as informações das vias identificadas e descartadas na análise
total <- nrow(table_information)
selected <- length(which(table_information$funn_pathways == 1))
exclude <-total - selected
slices <- c(exclude, selected)
lbls <- c(paste0("Excluded (", exclude,")"), paste0("Selected (", selected,")"))

## PT[Br] Retirando informações de 'global and overviews maps' [-g desativado] das vias selecionadas nas análises estatísticas
funn_indexes <- which(table_information$funn_pathways == 1)
global_indexes <- which(gsub(" ","\\1",as.character(table_information$subclass)) == "Globalandoverviewmaps")
if(S == TRUE){funn_indexes <- setdiff(funn_indexes,global_indexes)}
if(output_status == TRUE) { 
  write.csv(table_information, file = paste0(OUTPUT_PATH,sample,"-pathways-output.csv"), row.names = F) 
  jpeg(file = paste0(OUTPUT_PATH,sample,"-group-output.jpeg"))
  pie(slices, 
      labels = lbls, 
      main= paste0("Pie Chart of ", sample))
  dev.off()
}

## PT[Br] Obtendo informações dos grupos ortólogos envolvidos na rede
ko_indexes <- which(grepl("path:map",as.character(rownames(table_metrics_mgpn))) == FALSE)
ko_informations <- table_metrics_mgpn[ko_indexes,]
names <- rownames(ko_informations)

## PT[Br] Capturando as funções dos grupos ortólogos via KEGGREST e escrevendo as informações em um arquivo .csv
ko_functions <- keggKos(rownames(ko_informations))
ko_functions <- as.data.frame(ko_functions)

ko_functions <- as.character(ko_functions$functions[1:nrow(ko_informations)])
ko_informations <- cbind(ko_informations, ko_functions, names)
colnames(ko_informations)[which(colnames(ko_informations) == "names")] <- "kegg_ko"
ko_informations <- as.data.frame(ko_informations)

ord <- order(forNumeric(ko_informations$load_roume), decreasing = TRUE)
ko_informations <- ko_informations[ord,]


## PT[Br] Setando o tamanho dos vértices da rede de acordo com a categoria e a métrica de betweennes centrality
min <- 20; max <- 200

size_vertex <- forNumeric(table_information$load_roume[funn_indexes])

df2 <- as.character(gp_df2$name[which(gp_df2$type != 'KO families')])
t_nfo <- as.character(table_information$pathways[funn_indexes])
index <- NULL
for(i in 1:length(df2)){
  index[i] <- which(t_nfo == df2[i])
}

n_groups <- length(which(gp_df2$type == 'KO families'))
size_groups <- rep(min,n_groups)
paths <- minMax(size_vertex[index],min,max)
gp_df2$nodoSize <- c(size_groups,paths)

gp_df2$Class <- as.factor(c(rep('KO families',n_groups),as.character(table_information$class[funn_indexes])))
gp_df2$Subclass <- as.factor(c(rep('KO families',n_groups),as.character(table_information$subclass[funn_indexes])))

id <- seq(1:length(gp_df2$name))
label <- as.character(gp_df2$name)
group <- as.character(gp_df2$Class)
nodes <- cbind(id,label,group)
nodes <- as.data.frame(nodes)

k <- 1; from <- NULL; to <- NULL
for(i in 1:length(nodes$label)){
  v1 <- as.character(nodes$label[i])
  size <- length(neighbors(g2,v1))
  for(j in 1:size){
    v2 <- neighbors(g2,v1)[[j]]$name
    index_v2 <- which(v2 %in% nodes$label)
    from[k] <- i
    to[k] <- index_v2 
    k <- k +1
  }
}

edges <- cbind(from,to)
edges <- as.data.frame(edges)

CONFIG_PATH <- gsub("/src","/setup.conf",getwd())
config_file <- read.table(file = CONFIG_PATH, sep = "=", header = FALSE)
OUTPUT <- as.character(config_file$V2[3])
OUTPUT <- gsub(" ","\\1",OUTPUT)
size <- nchar(OUTPUT)
if(substr(OUTPUT,size,size) != "/"){ OUTPUT <- paste0(OUTPUT,"/")}

project <- toupper(PROJECT_NAME)
sample <- toupper(SAMPLE_NAME)

save.image(file = paste0(OUTPUT, "funn-", project,"-",sample,".RData"))
#load("/home/leandro/R/version/funn-mg-v.13/data/tempfiles/funn-CANGA-MG_34_EMMA.RData")
cat("\n\nInserting the data into the database...\n\n\n")
## PT[Br] Copiando as informações pro banco de dados Metagnômico (MONGO)
library(mongolite)
library(jsonlite)

dbname <- "funn"
sampledb <- "manager"
host <- gsub(" ","\\1",as.character(config_file$V2[1]))
port <- gsub(" ","\\1",as.character(config_file$V2[2]))

conect_mongo <- paste0("mongodb://", host, ":", port)

comand <- paste0("python ",
                  gsub("/src","/python_db/index_all.py",getwd()),
                  " -s ", sample,
                  " -p ", project)
system(comand, wait = TRUE)

collection <- paste0(project,"-",sample)

mongo = mongo(collection = collection, db = dbname, url = conect_mongo,
              verbose = FALSE, options = ssl_options())

sample_tool = mongo(collection = "samples", db = sampledb, url = conect_mongo,
                    verbose = FALSE, options = ssl_options())

updat <- sample_tool$count(query = paste0('{"$and":[{"project" : "',project,'",
                                                     "sample_name":"',sample,'"}]}'))

if(updat == 0){
  
  sample_tool$insert(paste0('{"sample_name":"',sample,'",
                              "project":"',project,'"}'))
  
}

rownames(ko_informations) <- NULL

for(i in 1:nrow(ko_informations)){
  ko <- as.character(ko_informations[i,"kegg_ko"])
  
  updat <- mongo$count(query = paste0('{"kegg_ko" :"', ko,'"}'))
  
  if(updat == 0){mongo$insert(ko_informations[i,])}
  
  p <- toJSON(names(as.factor(neighbors(g2,ko))), pretty=TRUE)
  toJSON(gp_df2)
  
  mongo$update(query = paste0('{"kegg_ko" :"', ko,'"}'),
               update = paste0('{"$set":{"paths": ',p,'}}'))
}

funn <- which(table_information$funn_pathways == 1)
pathway_information <- table_information[funn,]
cexc <- which(colnames(pathway_information) == "funn_pathways")
pathway_information <- pathway_information[,-cexc]

pathway_information$betweenness_centrality <- minMax(pathway_information$betweenness_centrality,0,100)

for(i in 1:nrow(pathway_information)){
  pathway <- as.character(pathway_information$pathways[i])
  
  updat <- mongo$count(query = paste0('{"pathways" :"', pathway,'"}'))
  
  if(updat == 0){mongo$insert(pathway_information[i,])}
  
  p <- toJSON(names(as.factor(neighbors(g2,pathway))), pretty=TRUE)
  
  mongo$update(query = paste0('{"pathways" :"', pathway,'"}'),
               update = paste0('{"$set":{"kos": ', p,'}}'))
  
}

sample_tool$update(query = paste0('{"$and":[{"project": "', project,'",
                                    "sample_name": "', sample,'"}]}'),
                   update = '{"$set":{"funn_tool": "funn-mg"}}')

if(PROTEOMIC == TRUE){
  
  comand <- paste0("python ",
                   gsub("/src","/python_db/proteomic-hit_mongo.py",getwd()),
                   " -i ", INPUT_PROT,
                   " -s ", SAMPLE_NAME,
                   " -p ", PROJECT_NAME)
  system(comand, wait = FALSE)
  
}

comand <- paste0("python ",
                 gsub("/src","/python_db/kaas_mongo.py",getwd()),
                 " -i ", INPUT_PATH,
                 " -s ", SAMPLE_NAME,
                 " -p ", PROJECT_NAME)
system(comand, wait = FALSE)

comand <- paste0("python ",
                 gsub("/src","/python_db/kaiju_mongo.py",getwd()),
                 " -i ", INPUT_TAXON,
                 " -s ", SAMPLE_NAME,
                 " -p ", PROJECT_NAME)
system(comand, wait = TRUE)


cat("\n\nExecution completed successfully.")
