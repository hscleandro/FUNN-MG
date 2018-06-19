library(RedeR)
library(igraph)
library(jsonlite)
library(mongolite)

CONFIG_PATH <- gsub("/src","/setup.conf",getwd())
config_file <- read.table(file = CONFIG_PATH, sep = "=", header = FALSE)

host <- gsub(" ","\\1",as.character(config_file$V2[1]))
port <- gsub(" ","\\1",as.character(config_file$V2[2]))

dbname <- "proteomics"
db_kaas <- "kaas"

project <- PROJECT_NAME
sample <- SAMPLE_NAME

conect_mongo <- paste0("mongodb://", host, ":", port)

collection <- paste0(project,"-",sample)

prot = mongo(collection = collection, db = dbname, url = conect_mongo,
              verbose = FALSE, options = ssl_options())
kaas = mongo(collection = collection, db = db_kaas, url = conect_mongo,
             verbose = FALSE, options = ssl_options())
#TYPE <- "pathway"

OUTPUT <- as.character(config_file$V2[3])
OUTPUT <- gsub(" ","\\1",OUTPUT)
size <- nchar(OUTPUT)
if(substr(OUTPUT,size,size) != "/"){ OUTPUT <- paste0(OUTPUT,"/")}

load(file = paste0(OUTPUT, "funn-", project,"-",sample,".RData"))

## PT[Br] Identificando o tipe de visualização escolhida pelo usuário
if(TYPE == 'class'){

  kegg_hier <- gp_df2$Class
  col <- "Class"
}else if(TYPE == 'subclass'){

  kegg_hier <- gp_df2$Subclass
  col <- "Subclass"
}else{
  kegg_hier <- gp_df2$type
  col <- "type"
}

if(PROT == TRUE){
  ko <- toJSON(gp_df2$name)
  prot_id <- prot$find(query = '{}',
                       fields = '{"id_seq":""}')
  id_seq <- toJSON(prot_id$id_seq)
  
  kaas_id <- kaas$find(query = paste0('{"$and":[{"id_seq":{"$in":', id_seq,'}, "kegg_ko":{"$ne":"nan"}}]}'),
                       field = '{"kegg_ko":""}')
  ko_prot <- kaas_id$kegg_ko
  ko_prot <- unique(ko_prot)
  
  index <- which(as.character(gp_df2$name) %in% ko_prot)
  gp_df2 <- as.matrix(gp_df2)
  gp_df2[index,col] <- "Proteomic"
  gp_df2 <- as.data.frame(gp_df2)
  kegg_hier <- gp_df2[,col]

}

## PT[Br] Selecionando as cores do vértices de acordo com as categorias
sizeColors <- length(unique(kegg_hier))
color <- rev(terrain.colors(sizeColors))
color_ko <- color[1]
ord <- unique(kegg_hier)[order(unique(kegg_hier), decreasing = T)]
n <- which(ord == 'KO families')
temp <- color[n]
color[n] <- color_ko
color[1] <-  temp
color <- rev(color)
if(PROT == TRUE){
  color_inst <- as.character(unique(kegg_hier))
  index <- order(color_inst, decreasing = FALSE) 
  color_inst <- color_inst[index]
  n <- which(color_inst == "Proteomic")
  color[n] <- "red"
}
if(TYPE == "pathways"){
  color_inst <- as.character(unique(kegg_hier))
  index <- order(color_inst, decreasing = FALSE) 
  color_inst <- color_inst[index]
  n <- which(color_inst == "Pathway")
  color[n] <- "#006600"
  
}
## PT[Br] setando parâmetros da visualização em rede

## PT[Br] Atribuindo características de nome, cor e tamanho ao objeto g1, considerando os dados de  gp_df
## Assigning name features, color and size to the object g1, considering the data  gp_df
g2 <- att.mapv(g2, dat= gp_df2, refcol=1)
g2 <- att.setv(g2, from="name", to ="nodeAlias")
g2 <- att.setv(g2, from=col, to="nodeColor", pal = 2, col= color)
g2 <- att.setv(g2, from="nodoSize", to="nodeSize", isrev=TRUE, 
             xlim=c(80,20,1))

#load('/home/leandro/Documents/teste_canga.RData')
## PT[Br] Visualizando a rede
rdp <- RedPort ()
calld(rdp)
resetd(rdp)
addGraph(rdp, g2, zoom=50.0, gscale=80, theme='tm1', layout.fruchterman.reingold(g2) )

addLegend.color(rdp, g2, position='bottomleft',vertical=TRUE, title="")
relax(rdp)

