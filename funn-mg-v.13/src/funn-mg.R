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

args <- commandArgs(trailingOnly = TRUE)

if("--help" %in% args){
  
  cat('\nFUNN-MG: FuNctioNal Network Analysis of MetaGenomics Data\n')
  cat('\nFUNN-MG should be executed in two sequential steps: (1) sequence analysis and (2) plot network.\n')
  cat('\n##SEQUENCE ANALYSIS\n')
  cat('\nHOW TO USE: (1) For analysis sequence do: Rscript funn-mg.R -s[required] -p[required] -f[required] -t[required] -o[optional] --pvalue[optional] --global[optional] --prot[optional]')
  cat('\n\nPARAMETER OPTIONS (1):\n')
  cat('\n -s: sample name (string)')
  cat('\n -p: project name (string)')
  cat('\n -f: Functional file of KAAS output (file path)')
  cat('\n -t: Taxonomic file of Kaiju output (file path)')
  cat("\n -o: output containg pathways identified in the sample (folder path)")
  cat("\n --pvalue: threshold of pvalue cut [default 0.05] (numeric 0..1)")
  cat("\n --global: consider Global and overview maps (not parameter) \n\n\n")
  cat('\n##PLOT NETWORK\n')
  cat('\nHOW TO USE: (2) For Network display do: Rscript funn-mg.R -s[required] -p[required] --display[required] --type[optional] --prot[optional]')
  cat('\n\nPARAMETER OPTIONS (2):\n')
  cat('\n -s: sample name (string)')
  cat('\n -p: project name (string)')
  cat('\n --display: Display network visualization (not parameter)')
  cat('\n --type: type of hierarchical visualization (class, subclass, pathways[default])')
  cat('\n --prot: fasta protein sequence file. (In --display mode does not need parameter)\n')
  
}else if("--display" %in% args){
  
  if("-f" %in% args){ cat("Erro: Parameter -f not requires the --display option.\n\n") }
  if("-t" %in% args){ stop("Erro: Parameter -t not requires the --display option.\n\n") }
  if("-o" %in% args){ stop("Erro: Parameter -o not requires the --display option.\n\n") }
  if("--pvalue" %in% args){ stop("Erro: Parameter --pvalue not requires the --display option.\n\n") }
  if("--global" %in% args){ stop("Erro: Parameter --global not requires the --display option.\n\n") }
  
  if("-s" %in% args){
    arg.sample <-  which(args == "-s") + 1
    SAMPLE_NAME <- as.character(args[arg.sample])
  }else{
    stop("Erro: Parameter -s required for execution.\n\n")
  }
  
  if("-p" %in% args){
    arg.project <-  which(args == "-p") + 1
    PROJECT_NAME <- as.character(args[arg.project])
  }else{
    stop("Erro: Parameter -p required for execution.\n\n")
  }
  
  if("--type" %in% args){
    arg.type <-  which(args == "--type") + 1
    TYPE <- as.character(args[arg.type])
  }else{
    TYPE <- "pathways"
  }
  
  if("--prot" %in% args){
    PROT <- TRUE
  }else{
    PROT <- FALSE
  }
  
  source(paste0(getwd(),"/display-grp.R"))
  
}else{
  
  system("clear")
  
  #if("--prot" %in% args){ stop("Erro: Parameter --prot requires the --display option.\n\n") }
  if("--type" %in% args){ stop("Erro: Parameter --type requires the --display option.\n\n") }
  
  if("-f" %in% args){
    arg.input <-  which(args == "-f")  + 1
    INPUT_PATH <- as.character(args[arg.input]) 
  }else{
    stop("Erro: Parameter -f required for execution.\n\n")

  }
  
  if("-t" %in% args){
    arg.taxon <-  which(args == "-t") + 1
    INPUT_TAXON <- as.character(args[arg.taxon])
  }else{
    stop("Erro: Parameter -t required for execution.\n\n")
  }
  
  if("-s" %in% args){
    arg.sample <-  which(args == "-s") + 1
    SAMPLE_NAME <- as.character(args[arg.sample])
  }else{
    stop("Erro: Parameter -s required for execution.\n\n")
  }
  
  if("-p" %in% args){
    arg.project <-  which(args == "-p") + 1
    PROJECT_NAME <- as.character(args[arg.project])
  }else{
    stop("Erro: Parameter -p required for execution.\n\n")
  }
  
  
  if("-o" %in% args){
    arg.output <-  which(args == "-o") + 1
    output_status <- TRUE
    OUTPUT_PATH <- as.character(args[arg.arg.output])
    OUTPUT_PATH <- gsub(" ","\\1",OUTPUT)
    size <- nchar(OUTPUT_PATH)
    if(substr(OUTPUT_PATH,size,size) != "/"){ OUTPUT_PATH <- paste0(OUTPUT_PATH,"/")}
  }
  else{ output_status <- FALSE }
  
  
  
  if("--prot" %in% args) { 
    arg.prot <-  which(args == "--prot") + 1
    INPUT_PROT <- as.character(args[arg.prot]) 
    PROTEOMIC <- TRUE
  }
  else{ PROTEOMIC <- FALSE }
  
  if("--pvalue" %in% args) {
    arg.pvalue <-  which(args == "--pvalue") 
    PVALUE <- as.numeric(as.character(args[arg.pvalue+1]))
  }
  else { PVALUE <- 0.05 }
  
  
  if("--global" %in% args) { 
    arg.global <-  which(args == "--global")
    S <- FALSE 
  }
  else { S <- TRUE }

  #view <- FALSE
  #INPUT_PATH = '/home/leandro/Data/BWWT/posdata/kaas/SRX612782_ghostx.kaas'
  #INPUT_TAXON = '/home/leandro/Data/BWWT/posdata/kaiju/SRX612782.kaiju.out.names'
  #SAMPLE_NAME = "SRX612782"
  #PROJECT_NAME = "BWWT"
  #PVALUE = 0.05
  #S <- FALSE
  
  DATASET <- strsplit(INPUT_PATH,"/")
  DATASET <- DATASET[[1]][length(DATASET[[1]])]
  #DATASET <- gsub('.kaas', "\\1", DATASET)
  cat("-------------------------------------------------------------------------\n")
  cat("******* FUNN-MG: FuNctioNal Network Analysis of MetaGenomics Data ******* \n")
  cat("Author: Leandro Correa - hscleandro@gmail.com\n")
  cat("        Ronnie Alves   - alvesrco@gmail.com \n")
  cat("Version: v1.13 - Copyright 2017 Leandro Corrêa\n")
  cat("-------------------------------------------------------------------------\n\n")
  
  
  ## Verifies if all the required packages are installed 
  source(paste0(getwd(),"/supports/packages.R"))
  ## library of functions that help during the execution
  source(paste0(getwd(),"/supports/bibfunction.R"))
  
  #tempfile <- paste0(gsub("/src","\\1",getwd()),"/data/temp.txt")

  cat("Input File:", DATASET,'\n')
  cat("P value:", PVALUE,"\n")
  cat("Type of analyse: KO Families")

  cat("\n\n\n")
  cat("Loading...\n\n\n")
  
  source(paste0(getwd(),"/pipeline-grp.R"))
  
  #source(paste0(getwd(),"/network_Metrics.R"))
  #source(paste0(getwd(),"/funn_Analysis.R"))
  
}





