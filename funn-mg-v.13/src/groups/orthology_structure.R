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
#INPUT_PATH <- '/home/leandro/Data/metagenomas/jcvi_samples/input/AM1/JCVI-AM1.kaas'
## PT[Br] Obtem a lista de grupos ortólogos no padrao KEGG
## Gets the list of orthology groups in KEGG standard
SAMPLE <- kaas_dataFrame(INPUT_PATH)

## PT[Br] Retira os grupos duplicados
## Remove duplicate groups
KO <- deduplicate(SAMPLE$KO)
KO <- gsub(" ","\\1",as.character(KO))

## PT[Br] Busca as vias metabólicas relacionadas a cada grupo
## Search pathways related to each group
temp <- keggPathways(KO)
rep <- which(substr(as.character(temp), 6, 7) == "ko")
pattern <- temp[-rep]

## PT[Br] Identifica vias e grupos
## Identify pathwyas and groups
paths <- as.character(pattern)
kos <- substr(as.character(names(pattern)),4,9)

## PT[Br] Organiza o nome das vias metabólicas para efetuar consulta no banco de dados KEGG
## Organizes the name of the metabolic pathways to make query in the KEGG database
group_path <- paths
names(group_path) <- substr(as.character(names(pattern)),4,9)

## PT[Br] Construindo a matriz KO grupos x pathways (1 onde existe relação grupo e via, 0 onde não existe)
## Building the KO groups x pathways matrix (1 where there is a relationship group and pathway, 0 where there isn't)
gp_MATRIX <- buildGPMatrix(group_path)


## PT[Br] Buscando a lista de todos os grupos relacionados aos pathways identificados na amostra e anotado no banco de dados KEEG
## Searching a list of all KO groups related to pathways identified in the sample and noted in the KEGG database
paths <- colnames(gp_MATRIX); kos <- rownames(gp_MATRIX)
kos_noted <- list(); sample <- NULL; noted <- NULL; total_kos <- NULL
for(i in 1:length(paths)){
  #newPath <- paste0("ko",split_str_by_index(paths[i],9)[2])
  #groups <- DB$group[which(as.character(DB$path) %in% newPath)]
  kos_noted[[i]] <- KeggLink("ko",paths[i])
  noted[i] <- length(kos_noted[[i]])
  total_kos <- c(total_kos, as.character(kos_noted[[i]]))
  sample[i] <- sum(gp_MATRIX[,i])
  #cat(i,"\n")
}

pathways <- colnames(gp_MATRIX)
pathway_information <- as.data.frame(cbind(pathways,sample,noted))
total_groups <- length(unique(total_kos))




