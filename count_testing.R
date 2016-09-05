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
# Rscript data_prep.R /path/counts_file /path/taxonomy_file /path/metadata.txt min_percentage_OTU_presence[0-1] /path/for/out_file/
# Rscript data_prep.R /path/16s.an.shared /path/16s.an.cons.taxonomy /path/metadata.txt 0.05 /path/for/out_file/
#
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
source(toolbox)
packages(c("metagenomeSeq"))

###### end ######

#* input *

f <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/results/dataF.rds' # commandArgs()[6] #
v <- "Salinity_InterstitialWater"#commandArgs()[7]
o <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/results/'#commandArgs()[8]

## Differential abundance testing ##
df <- df <- readRDS(f)
df <- cumNorm(df, p=0.5)
s <- normFactors(df)
mod <- model.matrix(~1 + pData(df)[[v]])
df.mod <- fitFeatureModel(df, mod)
head(MRcoefs(df.mod))

names(df.mod)









