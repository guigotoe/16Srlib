####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: March 2016
# Created: 29 September 2016
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript data_prep.R -h
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
#toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/age_lib/toolbox.R'
toolbox <- "/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse","ggplot2","outliers","devtools"))

## Options ##
#p = '/home/torres/ikmb_storage/Metagenome/16s/03.2016/'#
#r <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib_test/age'
p <- '/home/torres/ikmb_storage/Mangrove/ITS/pipits/single/F/ITS1/out_process/'#'/home/torres/ikmb_storage/Mangrove/16Sfa/mothur/'
r <- '/home/torres/Documents/Projects/Mangrove/ITS/results/'#'/home/torres/Documents/Projects/Mangrove/Results/16s/'

option_list <- list(
  make_option(c("-B","--biom"),action="store",type="character",default=paste(p,'otu_table.biom',sep=''),#'mgv16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.biom.bkup',sep=''),
              help="Path to input biom file"),
  make_option(c("-c","--counts"),action="store",type="character",default=paste(p,'otu_table.txt',sep=''),#paste(p,'mgv16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared',sep=''),
              help="Path to input counts file"),
  make_option(c("-x","--taxonomy"),action="store",type="character",default=NA,#paste(p,'mgv16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy.bkup',sep=''),
              help="Path to input taxonomy file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(p,'metadata',sep=''),
              help="Path to input metadata file"),
  make_option(c("-n","--out"),action="store",type="character",default='ITSmgv2017',
              help="Out files name"),
  make_option(c("-o","--outpath"),action="store",type="character",default=paste(r,'',sep=''),
              help="Path to output directory"),
  make_option(c("-t","--shared"),type="double",default=0.05,
              help="OTU presence; percentage of samples sharing each OTU. 0-1; default: %default"),
  make_option(c("-d","--depth"),type="double",default=10,
              help="Minimum depth count; default: %default counts per otu"),
  make_option(c("-l","--libt"),type="double",default=1.3,
              help="Library size threshold; times of sd threshold, default: %default - 0 means no threshold"),
  make_option(c("-r","--replicas"),type="logical",default=FALSE,
              help="Relicated libraries; default: %default"),
  make_option(c("-P","--PIPITS"),type="logical",default=TRUE,
              help="file from pipits; default: %default")
)
parser <- OptionParser(usage = "1. %prog -b path/to/biomfile -m path/to/metadata -o path/to/outdir [options] \n
                       2. %prog -c path/to/countfile -x path/to/taxonomy -m path/to/metadata -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if(!opt$PIPITS){
  ifelse(!is.na(opt$biom),using_biom <- TRUE,ifelse(!is.na(opt$counts),using_counts <- TRUE,using_counts <- FALSE))
  if (using_counts){
    if(is.na(opt$taxonomy)) {stop(sprintf("There is not taxonomy file specified"))}
    if(is.na(opt$metadata)) {stop(sprintf("There is not metadata file specified"))}
  }else if(using_biom==FALSE & is.na(opt$metadata)) {stop(sprintf("There is not biom nor counts nor metadata file specified"))
  }else if(using_biom & is.na(opt$metadata)) stop(sprintf("There is not metadata file specified"))
}else if (is.na(opt$counts) | is.na(opt$metadata)) stop(sprintf("Either counts or metadata files are missing"))


if(length(grep("/$",opt$outpath))==0) opt$outpath <- paste(opt$outpath,"/",sep="")

###### end ######

#* input *

source(toolbox)
message("Preparing the files...")
if(opt$PIPITS) {data <- pipts.otu(opt$counts,opt$metadata)
}else ifelse(using_biom, data <- mothur.biom(opt$biom,opt$metadata), data <- mothur.usingcounts(opt$counts,opt$metadata,opt$taxonomy))

message("Filtering...")

### Sample filtering by library size ### 
if(opt$libt != 0) data <- libsizefilter(data,opt$libt,opt$outpath)

## Using replicas to set filtering thresholds. 
if(opt$replicas) data <- replicas.analysis(data,opt$out)


if(as.numeric(opt$shared)==0){shared <- 1}else shared <- round(as.numeric(opt$shared)*NROW(pData(data)))

counts.f1.n <- MRcounts(data,norm=T,log=T)

featuresToKeep <- c() # OTUs satisfying presence and coverage thresholds
for (i in 1:NROW(counts.f1.n)){
  x <- counts.f1.n[i,]
  x <- x[!x %in% c(0)] # take off all 0 counts
  ## Identifying which OTUs satisfying only coverage threshold => shared = 1
  #if(length(x)>=1 & mean(x)>=otu.coverage) featuresToKeep.a <- c(featuresToKeep.a,rownames(counts.f1.n)[i])
  ## Identifying which OTUs satisfying presence and coverage thresholds
  if(length(x)>=shared & mean(x)>=log(opt$depth)) featuresToKeep <- c(featuresToKeep,rownames(counts.f1.n)[i])
}
data.f <- data[featuresToKeep,]

## geting raw information##
retained.info <- sum(fData(data.f)$Size)/sum(fData(data)$Size)
fData(data.f) <- droplevels(fData(data.f)) ## eliminate taxonomic levels without representation 
phylum <- length(levels(as.factor(fData(data.f)$Phylum))); class <- length(levels(as.factor(fData(data.f)$Class)))
order <- length(levels(as.factor(fData(data.f)$Order))); family <- length(levels(as.factor(fData(data.f)$Family)))
genus <- length(levels(as.factor(fData(data.f)$Genus)));otu <- length(levels(as.factor(fData(data.f)$OTU)))
pData(data.f)$libsize <- apply(MRcounts(data.f),2,sum)

## Normalizing
message("Normalizing...")
p.f <- cumNormStatFast(data.f) # Calculates the percentile for which to sum counts up to and scale by.
data.f <- cumNorm(data.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(data.f)
# saving files
message("Exporting files...")
saveRDS(data.f,file=paste(opt$outpath,opt$out,'_dataF','.rds',sep=''))
saveRDS(MRcounts(data.f),file=paste(opt$outpath,opt$out,'_otu_counts.rds',sep=''))
saveRDS(pData(data.f),file=paste(opt$outpath,opt$out,'_phenotype.rds',sep=''))
saveRDS(fData(data.f),file=paste(opt$outpath,opt$out,'_taxonomy.rds',sep=''))
output <- paste(paste("Filtred mode ",opt$out,":\n- An OTU was removed when its average read count was less than ",opt$depth," and when it was present in less than the ",100*as.numeric(opt$shared),"% of the samples\n",
                "- From ",NROW(MRcounts(data))," OTUs, ",NROW(MRcounts(data.f))," OTUs remained\n",
                "- Information (Reads) retained: ",round(100*retained.info,2),"%; Lost: ",round((100-100*retained.info),2),"%",sep=""),"\n",
                "- Taxonomic information:\n","   Phylum: ",phylum,"\tClass: ",class,"\tOrder: ",order,"\tFamily: ",family,"\tGenus: ",genus,"\tObserved OTU: ",NROW(MRcounts(data.f)),"\n",
                '- ',opt$out,"_dataF.rds -> Counts filtered and normalized - Cumulative-sum scaling normalization (Paulson et. al 2013)\n",
                "   \n** Data was successfully prepared! **\n",sep='')
message(output)
cat(output,file=paste(opt$outpath,'Report_data_prep.txt',sep=''),append=T)
