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

### Conection with KEGG and organization of information 
source(paste0(getwd(),"/groups/orthology_structure.R"))

## Functional enrichment analysis of metabolic pathways from the groups
source(paste0(getwd(),"/groups/orthology_enrichment.R"))

## For -f executions only the scripts bellow will run

## Build and extract important metrics of the network 
source(paste0(getwd(),"/network_Metrics.R"))
#save.image(paste0(temp_files,"/network_Metrics_",file))

## Visual representation and analysis of the FUNN-MG network about groups analysis
source(paste0(getwd(),"/funn_Analysis.R"))
