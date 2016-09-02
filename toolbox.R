####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: September 2016
# Created: 29 October 2015
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
#
# toolbox.R this script contains several functions
# commonly used by the 16Srlib pipeline for 16S analysis.
#
####################################################

packages <- function(requirements){
  has   <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=1:10)
    install.packages(requirements[!has])
  }
  lapply(requirements, require, character.only = TRUE)
}

numDFtranspose <- function(df){
  x <- t(df)
  colnames(x) <- x[1,]
  x <- x[-1,]
  z <- as.data.frame(apply(x, 2, as.numeric))
  return(z)
}

mothur.taxonomy <- function(taxonomy){
  packages(c("reshape2"))
  x <- colsplit(taxonomy$Taxonomy,pattern="\\([[:digit:]]*\\);",names = c('Kingdom','Phylum','Class','Order','Family','Genus'))
  x$Genus <- gsub("\\([[:digit:]]*\\);",'',x$Genus)
  x[x==""]  <- NA
  df <- cbind(taxonomy[,1:2],x)
  rownames(df) <- df[,1]
  return(df)
}
mothur.counts <- function(counts){
  counts$X <- NULL
  counts <- counts[,-c(1,3)]
  colnames(counts)[1] <- "id"
  df <- numDFtranspose(counts)
  rownames(df) <- rownames(df)
  return(df)
}
mothur.metadata <- function(metadata){
  rownames(metadata) <- metadata[,1]
  return(metadata)
}

###
###
# Geting sourcing paths
#initial.options <- commandArgs(trailingOnly = FALSE)
#file.arg.name <- "--file="
#script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
#script.basename <- dirname(script.name)
#other.name <- paste(sep="/", script.basename, "toolbox.R")
#print(paste("Sourcing",other.name,"from",script.name))
###

