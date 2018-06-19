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
if("gtable" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling gtable R package...\n")
  install.packages("gtable")
  #install.packages('../packages/Linux/gtable_0.2.0.tar.gz',repos=NULL, type="source")  
}
if("scales" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling scales R package...\n")
  install.packages("scales")
  #install.packages('../packages/Linux/scales_0.4.0.tar.gz',repos=NULL, type="source")  
}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling reshape2 R package...\n")
  install.packages("reshape2")
  #install.packages('../packages/Linux/reshape2_1.4.1.tar.gz',repos=NULL, type="source")  
}
if("XML" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling XML R package...\n")
  install.packages("XML")
  #install.packages('../packages/Linux/XML_3.98-1.4.tar.gz',repos=NULL, type="source")  
}
if("RCurl" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling RCurl R package...\n")
  install.packages("RCurl")
  #install.packages('../packages/Linux/RCurl_1.95-4.8.tar.gz',repos=NULL, type="source")  
}
if("qvalue" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling qvalue R package...\n")
  source("http://bioconductor.org/biocLite.R")
  biocLite("qvalue")
  #install.packages('../packages/Linux/Rcpp_0.12.8.tar.gz',repos=NULL, type="source") 
  #install.packages('../packages/Linux/colorspace_1.3-1.tar.gz',repos=NULL, type="source")
  #install.packages('../packages/Linux/plyr_1.8.4.tar.gz',repos=NULL, type="source")
  #install.packages('../packages/Linux/stringi_1.1.1.tar.gz',repos=NULL, type="source")
  #install.packages('../packages/Linux/qvalue_2.4.2.tar.gz',repos=NULL, type="source")  
}
if("mongolite" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling mongolite R package...\n")
  
  install.packages("mongolite")  
}

