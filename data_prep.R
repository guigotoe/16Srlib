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
# Rscript data_prep.R -c ../16Srlib_test/16S.otus.count -x ../16Srlib_test/16S.otus.taxonomy -m ../16Srlib_test/metadata -o ../16Srlib_test/results
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
#toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse"))

## Options ##
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'


option_list <- list(
  make_option(c("-c","--counts"),action="store",type="character",default=NA,#paste(p,'16S.otus.count',sep=''),
              help="Path to input counts file"),
  make_option(c("-x","--taxonomy"),action="store",type="character",default=NA,#paste(p,'16S.otus.taxonomy',sep=''),
              help="Path to input taxonomy file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=NA,#paste(p,'metadata',sep=''),
              help="Path to input metadata file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'results',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-t","--shared"),type="double",default=0.2,
              help="OTU presence; percentage of samples sharing each OTU. 0-1; default: %default"),
  make_option(c("-d","--depth"),type="double",default=10,
              help="Minimum depth count; default: %default counts per otu")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$counts)){stop(sprintf("There is not counts file specified"))
}else if(is.na(opt$taxonomy)){stop(sprintf("There is not taxonomy file specified"))
}else if(is.na(opt$metadata)){stop(sprintf("There is not metadata file specified"))}
if(length(grep("/$",opt$out))==0) out <- paste(opt$out,"/",sep="")


###### end ######

#* input *

message("Preparing the files...")
metadata <- mothur.metadata(read.table(opt$metadata,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
taxonomy <- mothur.taxonomy(read.table(opt$taxonomy,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- mothur.counts(read.table(opt$counts,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- subset(counts,taxonomy$Kingdom!="unclassified")
taxonomy <- subset(taxonomy,taxonomy$Kingdom!="unclassified")
#length(colnames(counts));length(union(colnames(counts),rownames(metadata)))
if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
  message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
  quit()
}

ord = match(colnames(counts),rownames(metadata))
metadata = metadata[ord,]
data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
message("Filtering...")
if(as.numeric(opt$shared)==0){shared <- 1}else{shared <- round(as.numeric(opt$shared)*NROW(pData(data)))}
data.f <- filterData(data,present=shared,depth=opt$depth)
retained.info <- sum(fData(data.f)$Size)/sum(fData(data)$Size)
message(paste(" - OTUs with less than ",opt$depth," counts and their presence in less than ",100*as.numeric(opt$shared),"% of the samples, were removed\n",
              " - From ",dim(MRcounts(data))[1]," OTUs, ",dim(MRcounts(data.f))[1]," OTUs remained\n",
              " - Information (Reads) retained: ",round(100*retained.info,2),"%; lost: ",round((100-100*retained.info),2),"%",sep=""))
## Normalizing
message("Normalizing...")
p.f <- cumNormStatFast(data.f) # Calculates the percentile for which to sum counts up to and scale by.
data.f <- cumNorm(data.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(data.f)
# saving files
message("Exporting files...")
saveRDS(data.f,file=paste(opt$out,'dataF.rds',sep=''))
message(" - dataF.rds -> Counts filtered and normalized - Cumulative-sum scaling normalization (Paulson et. al 2013)\n",
        " ** Data was successfully prepared! **")


