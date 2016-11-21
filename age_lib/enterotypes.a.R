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
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/age_lib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
packages(c("metagenomeSeq","reshape2","vegan","ggplot2","optparse","matrixStats","pheatmap","clusterSim"))

## Options ##
p <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib_test/age/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/age/'

option_list <- list(
  make_option(c("-i","--data"),action="store",type="character",default=paste(p,'dataFcp.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'',sep=''),#"./",#
             help="Path to output directory [default %default]")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not file specified"))}
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

#### Preparing the input data ####
df <- readRDS(opt$data)
pData(df) <- pData(df)[order(pData(df)$Age),]  # soring phenotype according Age
df <- df[,rownames(pData(df))]                 # soring samples according Age
pData(df)$group <- as.factor(unlist(apply(pData(df),1,function(x){
  if (as.numeric(x[["Age"]]) < 40) {return("G1")
  }else if (as.numeric(x[["Age"]]) > 39 & as.numeric(x[["Age"]]) < 60){ return("G2")
  }else if (as.numeric(x[["Age"]]) > 59 & as.numeric(x[["Age"]]) < 80){ return("G3")
  }else if (as.numeric(x[["Age"]]) > 79 & as.numeric(x[["Age"]]) < 90){ return("G4")
  }else if (as.numeric(x[["Age"]]) > 89) return("LLI")
})))

## Genus matrix
df.r <- newMRexperiment(returnAppropriateObj(df,norm=T,log=F),phenoData=AnnotatedDataFrame(pData(df)),featureData=AnnotatedDataFrame(fData(df)))
df.t <- aggTax(df.r,lvl="Genus")
df.f <- filterData(df.t,present=100)
df.g <- df.f[,which(pData(df.f)$group=="LLI")]
df.g <- df.f
## creating distance matrix
mat <- MRcounts(df.g,norm=T,log=F)
mat <- mat[-which(rownames(mat)=="unclassified"),]
data.dist=dist.JSD((mat))
data.dist=vegdist(t(mat))
nclusters=NULL
for (k in 1:8) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k=k)
    nclusters[k]=index.G1(t(mat),data.cluster_temp,d=data.dist,centrotypes="medoids")
  }
}
#pdf(paste(opt$out,"OptClusterNum_CHindex.pdf",sep=''),width=8, height=5)
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
#dev.off()
nc <- which(nclusters==max(nclusters[-c(1:2)]))
nc
data.cluster=pam.clustering(data.dist, k=nc)
head(data.cluster)

obs.pca=dudi.pca(data.frame(t(mat)),scannf=F,nf=10)
obs.bet=bca(obs.pca,fac=as.factor(data.cluster),scannf=F,nf=k-1)
plot(obs.bet)
s.class(obs.bet$ls,fac=as.factor(data.cluster),grid=F)


obs.coa=dudi.coa(data.frame(t(mat)),scannf=F)
bet <- bca(obs.coa,scannf=F)

## filtration ##
heatmapColColors <- brewer.pal(9, "Set1")[as.integer(factor(pData(df.f)$group))] 
heatmapCols <- colorRampPalette(c("#D9D9D9","#3794BF","#238B45","#D53E4F"))(4)
fontsize <- 6

df.f <- filterData(df.g,present=50)
matraw <- returnAppropriateObj(df.g,norm=T,log=T)
pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
         cluster_cols=F,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.g)[,("Age")])
pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
         cluster_cols=T,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.g)[,("Age")])





