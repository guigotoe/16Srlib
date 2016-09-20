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
# Rscript beta_div.R /path/dataF.rds association_variable co-variants/path/outfolder/
# Rscript beta_div.R ~/16Srlib_test/results/dataF.rds Salinity Limo,Arena ~/16Srlib_test/results/
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
p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'
packages(c("metagenomeSeq","reshape2"))

###### end ######

#* input *

f <- paste(p,'results/dataF.rds',sep='') #commandArgs()[6] # paste(p,'results/dataF.rds',sep='') #
vs <- 'Salinity' #commandArgs()[7]# 'Salinity' #
#vs <- unlist(strsplit(vs,','))
cf <- ''# commandArgs()[8]# # ,CT,NT,Ca,K,Mg,Na,CICE,Cu,S,P,Fe,Mn,Zn,B,Arcilla,Limo,Arena'
cf <- unlist(strsplit(cf,','))
#th <- 0.90 
o <- paste(p,'results/',sep='') #commandArgs()[9] # paste(p,'results/',sep='') #
## ##
df <- readRDS(f)
