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
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript apha_div.R /path/dataF.rds
#
#* requirements *#

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "toolbox.R")
source("/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R")
source(other.name)

packages(c("metagenomeSeq"))

###### end ######

#* input *

f <- '/home/torres/Documents/Projects/Metagenome/results/plotsMothur/09.2016/dataF.rds' # commandArgs()[6] #
df <- readRDS(f)

## diversity##
data(BCI)




