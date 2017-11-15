#!/usr/bin/env Rscript
####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: September 2016
# Created: September 2016
#
# This is written as part of 16S - Mangrove analysis, but could be
# splitted to serve different purposes.
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript beta_div.m.R -h
# Rscript beta_div.m.R -i ~/16Srlib_test/results/dataF.rds -o ~/16Srlib_test/results/ -e [options] 
#* requirements *#

get_script_path <- function() {
  cmdArgs = commandArgs(trailingOnly = FALSE)
  needle = "--file="
  match = grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    ls_vars = ls(sys.frames()[[1]])
    if ("fileName" %in% ls_vars) {
      # Source'd via RStudio
      return(normalizePath(sys.frames()[[1]]$fileName)) 
    } else {
      # Source'd via R console
      return(normalizePath(sys.frames()[[1]]$ofile))
    }
  }
}
script.basename <- dirname(get_script_path())
toolbox <- paste(sep="/", script.basename, "toolbox.R")
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
packages(c("metagenomeSeq","reshape2","vegan","ggplot2","optparse"))

## Options ##

p <- '/home/torres/Documents/Projects/Metagenome/2017_results/9_2017/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-i","--data"),action="store",type="character",default=paste(p,'dataFcp_l_0.2.rds',sep=''),
              help="Path to input rds file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'beta/',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-e","--exploratory"),action="store_true",default="NA",
              help="Perform exploratory analysis"),
  make_option(c("-m","--model"),action="store_true",default="NA",
              help="Build constrained model based on AIC selection criterion"),
  make_option(c("-am","--assess_model"),action="store",type="character",default=NA,#NULL,
              help="Model's terms assessed by permutation tests; for new model: vars,separated,by,comma"),
  make_option(c("-C","--constraints"),action="store",type="character",default=NULL,
              help="Set of constraints used by -coa: vars,separated,by,comma"),
  make_option(c("-coa","--constrained_analysis"),action="store_true",default="NA",
              help="Perform constrained ordination analysis using -C constraints"),
  make_option(c("-a","--factors"),action="store",type="character",default=NULL,
              help="Set of factors used by -b: vars,separated,by,comma"),
  make_option(c("-b","--beta"),action="store_true",default='NA',
              help="Perform beta diversty between -a variable classes"),
  make_option(c("-f","--filter"),action="store",type="double",default=0.3,
              help="Percentile as threshold of low abundant features ")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not file specified"))}

##
if(dir.exists(opt$out)){message('Out-folder already exist, files will be overwritten')
}else dir.create(opt$out,showWarnings=F)

#### Preparing the input data ####
data <- readRDS(opt$data)

packages('mvabund')
herb <- read.csv(file=paste(p,"Herbivore_specialisation.csv",sep=''),header=TRUE)
herbsp <- mvabund(herb[,5:11])
head(herbsp)
counts <- mvabund(t(MRcounts(data,norm=T)),check.names=F,row.names=colnames(MRcounts(data)))

meanvar.plot(counts,xlab='Mean',ylab='Variance')

pData(data)$age_group <- as.factor(unlist(apply(pData(data),1,function(x){
  if (as.numeric(x[["age"]]) <= 40) {return("G1")
  }else if (as.numeric(x[["age"]]) > 40 & as.numeric(x[["age"]]) <= 60){ return("G2")
  }else if (as.numeric(x[["age"]]) > 60 & as.numeric(x[["age"]]) <= 80){ return("G3")
  }else if (as.numeric(x[["age"]]) > 80) return("G4")
})))

plot(counts~pData(data)$age_group,cex.axis=0.8,cex=0.8)

mod1 <- manyglm(counts ~ pData(data)$age_group+pData(data)$gender, family="negative_binomial")
plot(mod1)
anovaTest

