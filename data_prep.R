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
#source("/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R")
#source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16Srlib/toolbox.R")
packages(c("metagenomeSeq"))

###### end ######

#* input *
c <- commandArgs()[6] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/16s.an.shared' # commandArgs()[6] #
t <- commandArgs()[7] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/16s.an.cons.taxonomy' # commandArgs()[7] #
m <- commandArgs()[8] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/metadata.txt' # commandArgs()[8] #
th <-commandArgs()[9] #0.05 #commandArgs()[9] # percentage threshold of OTU's presence across the samples.
d <- 10 # depth count threshold - by default.
o <- commandArgs()[10] #'/home/torres/Documents/Projects/Metagenome/results/plotsMothur/09.2016/' #


message("Preparing the files...")
metadata <- mothur.metadata(read.table(m,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
taxonomy <- mothur.taxonomy(read.table(t,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- mothur.counts(read.table(c,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
length(colnames(counts))
length(union(colnames(counts),rownames(metadata)))

if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
  message("**Error: Count-sample's names don't match with Metadata-sample's names!\n*Check the files and try again")
  quit()
}

ord = match(colnames(counts),rownames(metadata))
metadata = metadata[ord,]
data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
message("Filtering...")
data.f <- filterData(data,present=round(as.numeric(th)*NROW(pData(data))),depth=d)
retained.info <- sum(fData(data.f)$Size)/sum(fData(data)$Size)
message(paste(" - OTUs with less than ",d," counts and their presence in less than ",100*as.numeric(th),"% of the samples, were removed\n",
              " - From ",dim(MRcounts(data))[1]," OTUs, ",dim(MRcounts(data.f))[1]," OTUs remained\n",
              " - Information retained: ",round(100*retained.info,2),"%; lost: ",round((100-100*retained.info),2),"%",sep=""))
## Normalizing
message("Normalizing...")
p.f <- cumNormStatFast(data.f) # Calculates the percentile for which to sum counts up to and scale by.
data.f <- cumNorm(data.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(data.f)
# saving files
message("Exporting files...")
saveRDS(data.f,file=paste(o,'dataF.rds',sep=''))
message(" - dataF.rds -> Counts filtered and normalized - Cumulative-sum scaling normalization (Paulson et. al 2013)\n",
        " ** Data was successfully prepared! **")


