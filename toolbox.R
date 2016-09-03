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
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has],repos="https://cran.uni-muenster.de/")
  }
  lapply(requirements, require, character.only = TRUE)
}

numDFtranspose <- function(df){
  x <- t(df)
  colnames(x) <- x[1,]
  x <- x[-1,]
  rowid <- rownames(x)
  z <- as.data.frame(apply(x, 2, as.numeric))
  rownames(z) <- rowid
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
  return(df)
}
mothur.metadata <- function(metadata){
  if (is.na(table(duplicated(metadata[,1]))["TRUE"])){row.names(metadata) <- metadata[,1]
  } else {rownames(df) <- make.names(metadata[,1], unique=TRUE)}
  return(metadata)
}
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
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

