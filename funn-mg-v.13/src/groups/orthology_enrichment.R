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

library(qvalue)
library(igraph)

Enrichment <- NULL; p.value <- NULL
total_sample <- nrow(gp_MATRIX)
for(i in 1:nrow(pathway_information)){

  # A variável hitInSample recebe os grupos associados com a via 'i' na amostra
  # The variable hitInSample get the associated groups with the pathway 'i' in the sample
  hitInSample <- forNumeric(pathway_information$sample[i])

  # PT[Br] A variável total_hit_noted recebe o total de grupos encontrdos no banco de dados KEGG associados a via 'i'
  # The variable total_hit_noted get the total groups in the KEGG database associated with metabolic pathway 'i'
  total_hit_noted <- forNumeric(pathway_information$noted[i])
  
  # PT[Br] hitInPop recebe os grupos associados a via 'i' que estão anotados no banco de dados KEGG porém não são encontrados na amostra
  # hitInPop get the groups associated with metabolic pathway 'i' that are noted in the KEGG database but are not found in the sample
  hitInPop <- total_hit_noted - hitInSample
  
  # PT[Br] A variável failInSample recebe os grupos NÃO associados com a via 'i' na amostra
  # The variable hitInSample get the NOT associated groups with the pathway 'i' in the sample
  failInSample <- total_sample - hitInSample
  
  # PT[Br] A variável total_fail_noted recebe o total de grupos identificados no banco de dados KEGG NÃO associados a via 'i'
  # The variable total_fail_noted get the total groups in the KEGG database NOT associated with metabolic pathway 'i'
  total_fail_noted <- total_groups - total_hit_noted
  
  # PT[Br] failInPop recebe os groups NÃO associados a via 'i' que estão anotados no banco de dados KEGG porém não são encontrados na amostra
  # failInPop get the groups NOT associated with metabolic pathway 'i' that are noted in the KEGG database but are not found in the sample
  failInPop <- total_fail_noted - failInSample
  
  Enrichment[i] <- dhyper(hitInSample, total_hit_noted, total_fail_noted, total_sample)
  p.value[i] <- fisher.test(matrix(c(hitInSample, hitInPop, failInSample, failInPop), 2, 2), alternative='greater')$p.value;
  
}

cat("Total ko groups found in the sample:",total_sample,"\n")
cat("Total ko groups identified in the KEGG database:",total_groups,"\n\n\n")

# PT[Br] Calculando o erro do tipo II
##  Calculating the Type II error
qvalue <- qValue(p.value)
if(is.null(qvalue)) {q.value <- rep(0,length(p.value))
}else {q.value <- qvalue$qvalues}

# PT[Br] Salvando a informação obtidas nos testes estatísticos para o conjunto de vias metabólicas
## Saving the information obtained from the statistical tests for the set of metabolic pathways
PATH_DATABASE <- paste0(gsub("src","\\1",getwd()),"data/local_database/brite_hierarchies.csv")
brite_DB <- read.csv(file = PATH_DATABASE)

brite <- matrix(data = "NA", nrow = length(pathways), ncol = 3)
colnames(brite) <- c("class","subclass","Function")
rownames(brite) <- pathways; i <- 1

for (path in pathways) {
  aux <- which(brite_DB$ID == path)
  if(length(aux) == 0){
    brite[i,1] <- 'NA'
    brite[i,2] <- 'NA'
    brite[i,3] <- 'NA'
  }else{
    brite[i,1] <- as.character(brite_DB$CLASS[aux])
    brite[i,2] <- as.character(brite_DB$SUBCLASS[aux])
    brite[i,3] <- as.character(brite_DB$FUNCTION[aux])
  }
  i <- i +1
}

brite <- as.data.frame(brite)
table_information <-  cbind(pathways,as.character(brite$class),as.character(brite$subclass),as.character(brite[,3]))
namesTI <- c("pathways",colnames(brite))

table_information <- cbind(table_information,forNumeric(pathway_information$sample),forNumeric(pathway_information$noted),Enrichment,p.value, q.value)
colnames(table_information) <- c(namesTI,"genes_relation_sample","genes_relation_noted","coverage","p.value","q.value")

## PT[Br] Criando um data frame com as informações dos grupos e das vias da gp_Mtrix
## Creating a data frame with the information of groups and pathways of gp_Mtrix
gp_df <- information_matrix_ko(gp_MATRIX)
g1 <- graph.incidence(gp_MATRIX)

#rm(list=setdiff(ls(), c('g1','table_information','tempfile','gp_df',
#                        'file','temp_files','SAMPLE','DATASET')))

#save.image(paste0(temp_files,"/",file))
