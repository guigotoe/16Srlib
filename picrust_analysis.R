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
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/age_lib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse"))

## Options ##
#pa = '/home/torres/ikmb_storage/Metagenome/16s/2017/4_2017/'                   # AGING VERSION
#ra <- '/home/torres/Documents/Projects/Metagenome/2017_results/9_2017/'        # AGING VERSION
pa <- '/home/torres/ikmb_storage/Mangrove/16Sfa/2016/mothur/'                   # MANGROVE VERSION
ra <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/picrust/'          # MANGROVE VERSION

option_list <- list(
  make_option(c("-i","--input"),action="store",type="character",default=paste(pa,'greengenes_result/16s_2016_salinity.biom',sep=''),#'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared'
              help="Path to input biom input file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(pa,'metadata.txt',sep=''),
              help="Path to input metadata file"),
  make_option(c("-c","--contribution"),action="store",type="character",default=paste(ra,'ko_contributions.txt',sep=''),
              help="Path to contribution file"),
  make_option(c("-d","--diff"),action="store",type="character",default=paste(ra,'GenusDifAbund_conf__CountsTaxonomy.txt',sep=''),
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
}else if(is.na(opt$metadata)){stop(sprintf("There is not metadata file specified"))}
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

picrust_path <- '/home/torres/Bin/picrust-1.1.2/scripts/' # local
picrust_outpath <- paste(opt$out,'analisys/',sep='')

if(dir.exists(picrust_outpath)){message('picrustOut folder already exist, files will be overwritten')
}else dir.create(picrust_outpath,showWarnings=F)

###### end ######
packages(c('biomformat','ggplot2','dplyr','dendextend','stringdist'))
ggdata <- mothur.biom(opt$input,opt$metadata)
pData(ggdata) <- pData(ggdata)[,colSums(is.na(pData(ggdata)))==0] ## remove coluns with NAs
pData(ggdata)$ID_ref <- factor(as.character(pData(ggdata)$ID_ref),levels=c('low','med','high'))
p.f <- cumNormStatFast(ggdata) # Calculates the percentile for which to sum counts up to and scale by.
ggdata <- cumNorm(ggdata,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(ggdata)
head(fData(ggdata))
## data ##
tpdata <- as.data.frame(t(pData(ggdata)))
tpdata$header <- as.character(rownames(tpdata))
tpdata <- tpdata[,c(ncol(tpdata),1:(ncol(tpdata)-1))]
ko2lfse <- read.table(paste(ra,'ko_predicted_metagenomes2LEfSe.txt',sep=''),header=F,sep="\t",quote="")
ko2lfse <- ko2lfse[,1:(ncol(ko2lfse)-1)]
colnames(ko2lfse) <- colnames(tpdata)
ko2lfse <- rbind(tpdata[c('Punto','ID_ref','Fe','Arena','Salinity'),],ko2lfse)
ko2lfse[1:10,1:10]
write.table(ko2lfse,paste(ra,'ko_predicted_metagenomes2LEfSe_ready.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)

diff_genus <- read.table(opt$diff,header=T,sep="\t",quote="")$'Genus'
df_cont <- read.table(opt$contribution,header=T,sep="\t",quote="")
head(df_cont)
#max(df_cont$ContributionPercentOfSample)
length(levels(df_cont$Genus))
k <- df_cont[grep(as.character(diff_genus[4]),df_cont$Genus),]
k <- k[with(k,order(-ContributionPercentOfAllSamples)),]

#dendjw <- stringdistmatrix(strsplit(as.character(diff_genus[4]),'_')[[1]][1],df_cont$Genus,method='jw',useNames='strings')
#dendjw[1,c(151,162)]
#which(dendjw==max(dendjw))

#%>%as.dist %>% hclust(method='average') %>% as.dendrogram


ggplot(k,aes(x=Sample,y=ContributionPercentOfAllSamples,fill=Genus))+geom_boxplot()



