####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: November 2017
# Created: November 2017
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript picrust_analysis.R 
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
toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/age_lib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse"))

## Options ##
#pa = '/home/torres/ikmb_storage/Metagenome/16s/2017/4_2017/'                   # AGING VERSION
#ra <- '/home/torres/Documents/Projects/Metagenome/2017_results/9_2017/'        # AGING VERSION
#pa <- '/home/torres/ikmb_storage/Mangrove/16Sfa/2016/mothur_gg/'                   # MANGROVE VERSION
#ra <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/'          # MANGROVE VERSION
pa <- '/Users/guillermotorres/Documents/Proyectos/UAN_Manglares/29_08_ResultsMgv16s_v2/'   # MANGROVE VERSION Home
ra <- '/Users/guillermotorres/Documents/Proyectos/UAN_Manglares/29_08_ResultsMgv16s_v2/tax4fun/'   # MANGROVE VERSION Home

option_list <- list(
  make_option(c("-i","--input"),action="store",type="character",default=paste(pa,'dataFcp_l_0.1.rds',sep=''),#'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared'
              help="Path to input rds file"),
  make_option(c("-db","--database"),action="store",type="character",default=paste(pa,'SILVA123/',sep=''),
              help="Path to input metadata file"),
  make_option(c("-c","--contribution"),action="store",type="character",default=NULL,#paste(ra,'ko_contributions.txt',sep=''),
              help="Path to contribution file"),
  make_option(c("-k","--diff"),action="store",type="character",default=NULL,#paste(ra,'GenusDifAbund_conf__CountsTaxonomy.txt',sep=''),
              help="Path to differentially abundant Genus file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(ra,'',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-t","--shared"),type="double",default=0.2,
              help="OTU presence; percentage of samples sharing each OTU. 0-1; default: %default"),
  make_option(c("-d","--depth"),type="integer",default=1,
              help="Minimum depth count; default: %default counts per otu"),
  make_option(c("-L","--libraries"),type="integer",default=3,
              help="Minimum depth count per libraries [depth,j]: cutoff at 'depth' counts in at least j libraries."),
  make_option(c("-r","--replicates"),action="store_true",default=FALSE,
              help="Technical replications")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$input)){stop(sprintf("There is not counts biom-file specified"))
}else if(is.na(opt$database)){stop(sprintf("There is not data-base folder specified"))}
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

if(dir.exists(opt$out)){message('tax4fun folder already exist, files will be overwritten')
}else dir.create(opt$out,showWarnings=F)

packages('Tax4Fun')

df.r <- readRDS(opt$input)
counts <- MRcounts(df.r)
otus <- unlist(apply(fData(df.r)[,c('Kingdom','Phylum','Class','Order','Family','Genus')],1,function(x) paste(x,collapse=';')))
rownames(counts) <- as.character(paste(otus,';',sep=''))
df.x <- list(sampleNames=colnames(MRcounts(df.r)),otuTable=counts)

funp <- Tax4Fun(df.x,opt$database)
Tax4Fun(GN16SData,opt$database)
head(GN16SData$otuTable)
dim(counts)
length(colnames(MRcounts(df.r)))
