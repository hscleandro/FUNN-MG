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

# Manager R Package

# Execute: sudo apt-get install libxml2-dev
#          sudo apt-get install libcurl4-gnutls-dev
# JRE version >= 6
########################################
if("shiny" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling shiny R package...\n")
  install.packages("shiny")
  #install.packages('../packages/Linux/gtable_0.2.0.tar.gz',repos=NULL, type="source")  
}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling dplyr R package...\n")
  install.packages("dplyr")
  #install.packages('../packages/Linux/scales_0.4.0.tar.gz',repos=NULL, type="source")  
}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling ggplot2 R package...\n")
  install.packages("ggplot2")
  #install.packages('../packages/Linux/reshape2_1.4.1.tar.gz',repos=NULL, type="source")  
}
if("shinyjs" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling shinyjs R package...\n")
  install.packages("shinyjs")
  #install.packages('../packages/Linux/XML_3.98-1.4.tar.gz',repos=NULL, type="source")  
}
if("plotly" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling plotly R package...\n")
  install.packages("plotly")
  #install.packages('../packages/Linux/RCurl_1.95-4.8.tar.gz',repos=NULL, type="source")  
}
if("DT" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling DT R package...\n")
  install.packages("DT")
  }
if("mongolite" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling mongolite R package...\n")
  
  install.packages("mongolite")  
}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling jsonlite R package...\n")
  
  install.packages("jsonlite")  
}
if("devtools" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling devtools R package...\n")
  
  install.packages("devtools")  
}


