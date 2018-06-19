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

library(igraph)

#save.image(file = "/home/leandro/Data/metagenomas/MG_34_Emma/funn_indexes.RData")
#load("/home/leandro/Data/metagenomas/MG_34_Emma/funn_indexes.RData")

## PT[Br] Carregando as informações obtidas no banco de dados KEGG a partir da execução anterior
#source(paste0(getwd(),"/supports/bibfunction.R"))
#load(paste0(temp_files,"/",file))

#tempfile_table <- read.csv(tempfile, header = FALSE)

#OUTPUT = as.character(tempfile_table$V1[1])
#TYPE = as.character(tempfile_table$V1[2])
#PVALUE = forNumeric(tempfile_table$V1[3])
#S = as.logical(tempfile_table$V1[4])
#view = as.logical(tempfile_table$V1[5])

#system(paste0('rm ',tempfile))
table_information <- as.data.frame(table_information)


## PT[Br] Obtendo as vias metabólicas com enriquecimento inferior ao pvalue informao pelo usuário 
less_than_pvalue <- which(forNumeric(table_information$p.value) < PVALUE)
if(nrow(table_information) > 100){
  less_than_qvalue <- which(forNumeric(table_information$q.value) < PVALUE)
  less_than_pvalue <- union(less_than_pvalue, less_than_qvalue)
}

funn_pathways <- rep(0,nrow(table_information))
funn_pathways[less_than_pvalue] <- 1

g2 <- g1
paths_discarded <- as.character(table_information$pathways[which(funn_pathways == 0)])

## PT[Br] Removendo 'Global views maps' [opção -g desativada]
if(S == TRUE){
  global_overview <- as.character(table_information$pathways[which(gsub(" ","\\1",as.character(table_information$subclass)) == "Globalandoverviewmaps")])
  paths_discarded <- union(paths_discarded, global_overview)
}

## PT[Br] Criando uma nova rede a partir das novas vias selecionadas
g2 <- delete.vertices(g2,paths_discarded)
groups_discarded <- as.character(names(which(degree(g2) == 0)))
g2 <- delete.vertices(g2,groups_discarded)
vertex_discarded <- c(paths_discarded, groups_discarded)

vertex_not_discarded <- which((as.character(gp_df$name) %in% vertex_discarded) == FALSE)
gp_df2 <- gp_df[vertex_not_discarded,]
gp_df2$degree <- forNumeric(degree(g2))

colnames(gp_df2) <- c("name","degree","type","nodoSize")
gp_df2 <- as.data.frame(gp_df2)

paths_selected <- as.character(gp_df2$name[which(gp_df2$type == 'Pathway')])
index <- which(as.character(table_information$pathways) %in% paths_selected)

class <- as.character(table_information$class[index])
subclass <- as.character(table_information$subclass[index])

## PT[Br] Extraindo medidas de grau de conectividade (gene-pathway Network)
## Extracting degree of connectivity (gene-pathway Network)
indegree_g2 <- degree(g2,v=V(g2),mode="in")
indegree_g2_deg_dist <- degree.distribution(g1,cumulative=T,mode="in")
indegree_g2_pwl <- power.law.fit(indegree_g2)

## PT[Br] Extraindo medida de centralidade de intermediação (gene-pathway Network)
## Extracting betweenness centrality (gene-pathway Network)
bet_centrality_g2 <- betweenness(g2, v=V(g2), directed = FALSE, 
                                 nobigint = TRUE, normalized = FALSE)

######################################_Test_###########################################
#teste <- read.csv(file = "/home/leandro/Data/outros/teste.csv", header = F)
#ns <- c("m1","m2","m3","m4","e1","e2","e3","e4")
#colnames(teste) <- ns
#rownames(teste) <- ns
#g3 <- graph.adjacency(as.matrix(teste))

#bet_centrality_g3 <- betweenness(g3, v=V(g3), directed = FALSE, 
                                 #nobigint = TRUE, normalized = FALSE)

#indegree_g3 <- degree(g3,v=V(g3),mode="in")
#######################################################################################

total_bet <- sum(bet_centrality_g2)

## PT[Br] Calculando o 'load_score' a partir das informações de betweenness e vertex degree segundo Batista D et. all (2015)
total_edg <- length(E(g2))/2; roume <- NULL; rahman <- NULL; heads <- NULL; choke <- NULL
for(i in 1:length(bet_centrality_g2)){
  name <- names(bet_centrality_g2)[i]
  heads[i] <- name
  roume[i] <- (bet_centrality_g2[i]/total_bet) / (indegree_g2[i]/total_edg)
  rahman[i] <- log((bet_centrality_g2[i]/indegree_g2[i]) / (total_bet/total_edg))
  choke[i] <- bet_centrality_g2[i]
}
names(roume) <- heads
names(rahman) <- heads
names(choke) <- heads

## PT[Br] Craindo uma tabela com os resultados das métricas extraídas da rede para as vias metabólicas identificadas
## Creating a table with the results of the metrics extracted from network for metabolic pathways identified
table_metrics_mgpn <- cbind(bet_centrality_g2, indegree_g2, roume, rahman, choke)
colummNames <- c("betweenness_centrality","degree", "load_roume", "load_rahman", "choke_point")
colnames(table_metrics_mgpn) <- colummNames
rownames(table_metrics_mgpn) <- as.character(gp_df2$name)


colname <- c(colnames(table_information),colummNames)
index_paths <- which(grepl("path:map",as.character(rownames(table_metrics_mgpn))) == TRUE)

selected_paths <- which(as.character(table_information$pathways) %in% rownames(table_metrics_mgpn)[index_paths])

betweenness_centrality <- rep(0,nrow(table_information))
betweenness_centrality[selected_paths] <- round(forNumeric(table_metrics_mgpn[index_paths,"betweenness_centrality"]),2)


index_paths <- which(grepl("path:map",as.character(names(roume))) == TRUE)
load_roume <- rep(NA,nrow(table_information))
load_rahman <- rep(NA,nrow(table_information))
choke_point <- rep(NA,nrow(table_information))

load_roume[selected_paths] <- roume[index_paths]
load_rahman[selected_paths] <- rahman[index_paths]
choke_point[selected_paths] <- choke[index_paths]
#load_score <- minMax(load_score,0,100)
  
table_information <- cbind(table_information, betweenness_centrality, funn_pathways, load_roume, load_rahman, choke_point)

## PT[Br] Removendo 'Global views maps' [opção -g desativada]
if(S == TRUE){
  lines_GlobalOver <- which(gsub(" ","\\1",as.character(table_information$subclass)) == "Globalandoverviewmaps")
  table_information$funn_pathways[lines_GlobalOver] <- 0
}



