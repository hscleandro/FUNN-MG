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

# Library that supports this tool 
########################################

# PT[Br] Consulta as classes identificadas a cada via metabólica encontrada na amostra
# Consult identified classes in each metabolic pathway found in the sample

brite_Hierarchies <- function(pathways){
  PATH_DATABASE <- paste0(gsub("src","\\1",getwd()),"data/local_database/brite_hierarchies.csv")
  brite_DB <- read.csv(file = PATH_DATABASE)
  
  brite <- matrix(data = "NA", nrow = length(pathways), ncol = 3)
  colnames(brite) <- c("class","subclass","function")
  rownames(brite) <- pathways; i <- 1
  #brite <- as.data.frame(brite)
  for (path in pathways) {
    aux <- which(brite_DB$ID == path)
    brite[i,1] <- as.character(brite_DB$CLASS[aux])
    brite[i,2] <- as.character(brite_DB$SUBCLASS[aux])
    brite[i,3] <- as.character(brite_DB$FUNCTION[aux])
    i <- i +1
  }
  print(brite)
  cat("\n\n\n")
  #brite <- as.data.frame(brite)
  return(brite)
}

# PT[Br] Constrói a matriz gene-pathway
# BUild gene-pathway matrix
buildGPMatrix <- function(gene_path){
  ## PT[Br] Obtendo o numero de linhas e colunas da matriz segundo o numero de genes e pathways da amostra
  ## Getting the number of rows and columns of the matrix according to the number of genes and pathways of the sample
  listpathway <- names(table(gene_path))
  listGenes <- names(table(names(gene_path)))
  sizeCol <- length(listpathway); sizeRow <- length(listGenes)
  #procarioto_table <- kegg_information(LOCAL_KEGG_PATH); nOrg <- nrow(table_org)
  
  ## PT[Br] Iniciando a matriz e atribuindo nomes a linhas (genes) e colunas (pathways)
  ## Starting the array and assigning names to rows (genes) and columns (pathways)
  gp_Matrix <- matrix(data = 0, nrow = sizeRow, ncol = sizeCol)
  colnames(gp_Matrix) <- listpathway
  rownames(gp_Matrix) <- gsub(".*:(\\w+)", "\\1", listGenes)
  
  ## PT[Br] verifica para cada pathway do vetor listpathway, quais genes estão relacionados, cada verificação verdadeira
  # é marcado com '1' na gp_Mathttp://rest.kegg.jp/link/ko/M00010rix, onde as colunas são os pathways e as linhas os genes encontrados na amostra.
  
  ## checks for each pathway of listpathway vector, which genes are related, each real check is marked with '1'
  # in gp_Matrix, where the columns are the rows pathways and the genes found in the sample.
  gene_path_COPY <- gene_path
  for(i in 1:length(listpathway)){
    targetGenes <- names(which(gene_path_COPY == listpathway[i]))
    n <- which(gene_path_COPY == listpathway[i]); names(n) <- NULL
    targetGenes <- gsub(".*:(\\w+)", "\\1", targetGenes)
    
    index <- which(rownames(gp_Matrix) %in% targetGenes)
    gp_Matrix[index,i] <- 1
    
    gene_path_COPY <- gene_path_COPY[-n]
  }
  rownames(gp_Matrix) <- gsub(".*:(\\w+)", "\\1", listGenes)
  return(gp_Matrix)
}

# PT[Br] Escolhe cores em elementos de um vetor
# Choose color elements of a vector
colors <- function(teste){
  idColor <- terrain.colors(length(unique(teste)))
  names(idColor) <- unique(teste)
  color1 <- teste
  for(i in 1:length(idColor)){
    n <- which(names(idColor)[i] == teste)
    color1[n] <- idColor[i]
  }
  return(color1)
}


# PT[Br] Remove os genes duplicados da amostra
# Remove deduplicate genes of sample
deduplicate <- function(target){
  torfSet <- duplicated(target)
  setGenes <- target
  t <- which(torfSet == TRUE)
  if(length(t) > 0)
    setGenes <- target[-t]
  null <- which(setGenes == "NULL")
  if(length(null) > 0)
    setGenes <- setGenes[-null]
  return(as.character(setGenes))
}

## PT[Br]  Transforma o tipo factor em número
## Turn factor to numeric type
forNumeric <- function(tfactor){
  return(as.numeric(as.character(tfactor)))
}

# PT[Br] Retonra um data frame com as informações de genes e pathways após a fase de enriquecimento
# Returns a data frame with the information about genes and pathways after enrichment phase
information_matrix <- function(gp_Matrix, table_org, listGenes){
  
  ngenes <- nrow(gp_Matrix); npaths <- ncol(gp_Matrix); size <- ngenes + npaths
  gene_name <- rownames(gp_Matrix); path_name <- colnames(gp_Matrix); 
  table_name <-c(gene_name, path_name)
  
  # PT[Br] variável que armazenará o grau dos nodos
  # variable that stores the degree of the nodes
  degree <- rep(0, size)
  
  # Obtendo o id no Kegg para os organismos
  # Getting the id in Kegg for organisms
  label_gene <- gsub('[:].*', "\\1",listGenes)
  label_path <- rep("path",npaths)
  label <- c(label_gene, label_path)
  
  # PT[Br] Obtendo os nomes das espécies
  # Getting the species names
  name_org <- rep("Pathway",size) 
  for(i in 1:length(label_gene)){
    index <- which(label[i] == table_org[,"Category"]) 
    name_org[i] <- table_org[index,"Organisms"]
  }
  
  nodoSize <- rep(0,size)  
  
  # Escrevendo as informações obtidas em uma tabela
  # Writing the information obtained in a data frame
  gp_df = data.frame(table_name, label, degree, name_org, nodoSize)
  colnames(gp_df) <- c("name","idKegg","degree","nameOrg","nodoSize")
  return(gp_df)
}

# PT[Br] Retorna um data frame com as informações de genes e pathways após a fase de enriquecimento (análise de grupos ortólogos)
# Returns a data frame with the information about genes and pathways after enrichment phase (groups orthologys analyse)
information_matrix_ko <- function(gp_Matrix){
  ngroups <- nrow(gp_Matrix); npaths <- ncol(gp_Matrix); size <- ngroups + npaths
  group_name <- rownames(gp_Matrix); path_name <- colnames(gp_Matrix); 
  
  name <-c(group_name,path_name)
  
  degree <- rep(0,size)
  
  ko <- rep("KO families",ngroups)
  path <- rep("Pathway",npaths)
  nameOrg <- c(ko,path)
  
  nodoSize <- rep(0, size)
  
  gp_df = data.frame(name, degree, nameOrg, nodoSize)
  return(gp_df)
}

# PT[Br] Converte o output da ferramenta KAAS para um dataframe
kaas_dataFrame <- function(INPUT_PATH){
  SAMPLE <-read.csv(INPUT_PATH, header = FALSE, sep = "\t")
  if(ncol(SAMPLE) != 2){
    head <- NULL; ko <- NULL; k <- 1
    for(i in 1:nrow(SAMPLE)){
      temp <- as.character(SAMPLE$V1[i])
      size_head <- nchar(temp)
      key_char <- substr(temp,1,2)
      
      if((size_head <= 8)&&(key_char == "K0")){
        ko[k] <- temp
        k <- k +1
      }
      else{
        head[k] <- temp
        
      }
    }
    SAMPLE <- cbind(head[1:length(ko)], ko)
    SAMPLE <- as.data.frame(SAMPLE)
  }
  
  colnames(SAMPLE) <- c("ID","KO")
  
  return(SAMPLE)
}

# PT[Br] Retonra uma funçao a partir dos identificadores das vias metabolicas anotadas no KEGG
# Return the function from the recorded metabolic pathways identifiers in the KEGG database
kegg_fuctions <- function(pathways){
  pathways <- as.character(pathways)
  list <- NULL
  size <- length(pathways); 
  if(size < 100)
    limit.max <- size
  else limit.max <- 100
  limit.min <- 1
  while(length(list) < length(pathways)){

    list <- c(list, KeggList(pathways[limit.min:limit.max]))

    size <- length(pathways) - length(list)
    limit.min <- limit.max + 1
    if(size > 100){
      limit.max <- limit.max + 100
    }
    else {limit.max <- limit.max + size}
  }
  return(list)
  
}

# PT[Br] Retonra os grupos ortólogos identificadores a partir dos módulos consultados
kegg_groups <- function(modules){
  KEGG_URL <- "http://rest.kegg.jp/link/ko/"
  if(length(modules) == 1) 
  {
    SEARCH <- modules
  }else
  {
    SEARCH <-  paste0(modules, collapse = "+")
  
  }
  
  URL_SEARCH <- paste0(KEGG_URL, SEARCH)
  
  kos <- NULL; k <- 1
  for (i in readLines(URL_SEARCH)){
    kos[k] <- strsplit(i,"\t")[[1]][2]
    k <- k +1
  }
  
  kos <- gsub('.*[:]', "\\1", kos)
  
  return(kos)
}

# PT[Br] Retorna a tabela com todas as infromações sobre os procariotos do banco de dados KEGG
# Retuns the table with all informations about prokaryotes of KEGG database
kegg_information <- function(PATH){
  procarioto_table <-read.csv(PATH, header=T, sep=",")
  procarioto_table <- as.matrix(procarioto_table)
  return(procarioto_table)
}

# PT[Br] A mesma funçao da funçao keggLink do pacote KEGGREST
# The same function of keggLink function of KEGGREST R package
KeggLink <- function(target, source){
  KEGG_URL <- "http://rest.kegg.jp/link/pathway/"
  if(target == "ko"){
    KEGG_URL <- "http://rest.kegg.jp/link/ko/"
  }
  genes <- paste0(source, collapse = "+")
  URL_SEARCH <- paste0(KEGG_URL, genes)
  if(length(source) == 1){
    if(source == "pathway"){
      KEGG_URL <- "http://rest.kegg.jp/link/pathway/"
      URL_SEARCH <- paste0(KEGG_URL, target)
    }
  }
  gene <- NULL; path <- NULL; query <- NULL; k <- 1
  for (i in readLines(URL_SEARCH)){
    #cat(i,"\n")
    query <- strsplit(i,split = "\t")[[1]]
    gene[k] <- query[1]
    path[k] <- query[2]
    k <- k +1
  }
  names(path) <- gene
  
  return(path)
}

# PT[Br] Retonra as funções a nível de proteína a partir da consulta de grupos ortólogos
keggKos <- function(vector){
  max <- 100; size <- length(vector); kos <- NULL; functions <- NULL
  if(size > max){
    i <- 1; k <- max; p <- 1
    step <- size %/% max +1
    while(p <= step){
      if(p == step)
        functions <- KeggList(vector[i:size])
      else if(p < step)
        functions <- KeggList(vector[i:(i+k)])
      functions <- as.data.frame(functions)
      kos <- rbind(kos, functions)
      i <- i+k; p <- p +1; functions <- NULL
    }
  }
  else{
    kos <- KeggList(vector[1:size])
    kos <- as.data.frame(kos)
  }
  return(kos)
}

# PT[Br] A mesma funçao da funçao keggList do pacote KEGGREST
# The same function of keggList function of KEGGREST R package
KeggList <- function(source){
  KEGG_URL <- "http://rest.kegg.jp/list/"
  pathways <- paste0(source, collapse = "+")
  URL_SEARCH <- paste0(KEGG_URL, pathways)
  
  
  info <- NULL; path <- NULL; query <- NULL; k <- 1
  for (i in readLines(URL_SEARCH)){
    query <- strsplit(i,split = "\t")[[1]]
    path[k] <- query[1]
    info[k] <- query[2]
    k <- k +1
  }
  names(info) <- path
  
  return(info)
}

# PT[Br] Recebe um vetor de genes e retonra todos os pathways relacionados anotados no banco de dados kegg
# Receives an array of genes and returns all related pathways noted in KEGG database
keggPathways <- function(vector){
  max <- 100; size <- length(vector); pathways <- NULL; temp <- NULL
  if(size > max){
    i <- 1; k <- max; p <- 1
    step <- size %/% max +1
    while(p <= step){
      if(p == step)
        temp <- KeggLink("pathway", vector[i:size])
      else if(p < step)
        temp <- KeggLink("pathway", vector[i:(i+k)])
      pathways <- c(pathways, temp)
      i <- i+k; p <- p +1; temp <- NULL
    }
  }
  else{
    pathways <- KeggLink("pathway", vector[1:size])
  }
  return(pathways)
}

# PT[Br] Retorna a lista da quantidade de pathways por espécie encontrada na amostra
# Returns the list of the amount of pathways for species found in the sample
listPathway_perSpecies <- function(gene_path, table_org, nc){
  
  pathways <- names(table(gene_path))
  pathway_list_aux <- gene_path; 
  names(pathway_list_aux) <- gsub("(\\w+):.*", "\\1", names(gene_path))
  
  pathway_list_sample <- matrix(data = 0, ncol = nc, nrow = length(pathways))
  procarioto_table <- kegg_information(LOCAL_KEGG_PATH)
  
  for(i in 1:length(pathways)){
    index <- which(pathway_list_aux %in% pathways[i])
    lin <- table(names(pathway_list_aux[index]))
    n <- which(table_org[,"Category"] %in% names(lin))
    pathway_list_sample[i,n] <- lin
  }
  rownames(pathway_list_sample) <- pathways; colnames(pathway_list_sample) <- table_org[,"Category"]
  return(pathway_list_sample)
}

# PT[Br] Padronização de um vetor dado um intervalo mínimo e máximo
# Standardization of a vector given a minimum and maximum range
minMax <- function(vet, min=0, max=1){
  v_n <- NULL
  for(i in 1:length(vet)){
    v_n[i] <- min + (vet[i] - min(vet))/(max(vet) - min(vet)) *   (max-min)
    v_n[i] <- round(v_n[i], 3)
  }
  
  return(v_n)
}

# PT[Br] O usuário deve pressionar qualquer tecla pra prosseguir a execução
# The user must press any key to continue the execution
mywait <- function() {
  tt <- tktoplevel()
  tkpack( tkbutton(tt, text='Continue', command=function()tkdestroy(tt)),
          side='bottom')
  tkbind(tt,'<Key>', function()tkdestroy(tt) )
  
  tkwait.window(tt)
}

# PT[Br] Calcula a media de conectividade entre os vizinhos de cada nodo
# Calculates the average connectivity between neighbors of each node
neighborhoodConnectivity <- function(g2){
  neighborhood_gg<-array()
  #Para cada nó calcular the neighborhood connectivity
  for (i in 1: vcount(g2)){
    neighbors<-neighbors(g2,i)
    if(length(neighbors) > 0){
      neighborhoodConnectivity <- 0
      for (j in 1: length(neighbors)){
        neighborhoodConnectivity <- as.numeric(degree(g2,neighbors[j])) + neighborhoodConnectivity
      }
      neighborhood_gg[i] <- neighborhoodConnectivity/length(neighbors)
    }
    else{
      neighborhood_gg[i] <- 0
    }
  }
  return(neighborhood_gg)
}

# PT[Br] Estima o tamanho dos vétices da rede
# Estimates the size of the network vertices
nodoSize <- function(df){
  paths <- df[which((df$nameOrg == "Pathway") == TRUE),]
  ngenes <- nrow(df) - nrow(paths)
  minDegree <- min(paths$degree)
  nodSize <- c(rep(minDegree,ngenes),paths$degree)
  
  return(minMax(nodSize,20,round(4*(var(nodSize)))))
}

# PT[Br] Identifica quantas deleçoes são necessárias para cada via metabólica ficar sem conexões na rede gene-pathway
# Identifies how many deletions are required for each pathway runs out connections on the pathway-genes network 
nullTest <- function(g3, genes, paths,dels){
  #table_deletion <- nullTest(g2, genes, paths)
  #g3 <- g2
  table_deletion <- matrix(data = 0, ncol = 2, nrow = length(paths))
  table_deletion[,1] <- paths; 
  size_genes <- length(genes)
  #dels <- 100
  for(i in 1:as.numeric(size_genes)){
    #cat("1.\n")
    if(length(genes) > dels)
      {del <- sample(genes, dels)} #; cat("2\n")}
    if(length(genes) <= dels)
    {  del <- genes #; cat("3.\n")
    }
    #cat(i,":",del, "\n")
    #browser()
    n <- which(genes %in% del)
    genes <- genes[-n]
    g3 <- delete.vertices(g3, del)
    #cat("4.\n")
    for(j in 1:length(paths)){
      #cat("5.\n")
      size_path <- length(neighbors(g3,paths[j]))
      #cat("size_path: ", size_path, "\n")
      if((size_path == 0)&&(table_deletion[j,2] == 0))
        table_deletion[j,2] <- i #;cat("5.\n")
      #cat("6.\n")
    }
    #at("7.\n")
    #cat("genes",length(genes),"\n")
    #browser()
    if(length(genes) == 0){
      colnames(table_deletion) <- c("pathway name", "number of deletions")
      return(table_deletion)
    }
    
  }
  
}

# PT[Br] Função detecta outliers em um conjunto de dados e retorna os índices dos elementos
# Function detect outliers in a set of data and returns the indices of the elements
outliers <- function(data){
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  extreme.threshold.upper = (iqr * 3.5) + upperq
  index <- which(data > extreme.threshold.upper)
  return(index)
}

outliers_min <- function(data){
  #data <- prob
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  extreme.threshold.lower = lowerq - (iqr * 1.5)
  index <- which(data < extreme.threshold.lower)
  return(index)
}

qValue <- function(pval){
  return(tryCatch(qvalue(pval, lambda=0.5, robust=TRUE), error=function(e) NULL))
}



#===================================================================================================