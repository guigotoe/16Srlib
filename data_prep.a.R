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
toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse","ggplot2"))

## Options ##
#p = '/home/torres/ikmb_storage/Metagenome/16s/03.2016/'#
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-c","--counts"),action="store",type="character",default=paste(p,'16s.an.shared',sep=''),
              help="Path to input counts file"),
  make_option(c("-x","--taxonomy"),action="store",type="character",default=paste(p,'16s.an.cons.taxonomy',sep=''),
              help="Path to input taxonomy file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(p,'design.txt',sep=''),
              help="Path to input metadata file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'age',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-t","--shared"),type="double",default=0.05,
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
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

###### end ######

#* input *

message("Preparing the files...")
metadata <- mothur.metadata(read.table(opt$metadata,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
taxonomy <- mothur.taxonomy(read.table(opt$taxonomy,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- mothur.counts(read.table(opt$counts,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- subset(counts,taxonomy$Kingdom!="unclassified")
taxonomy <- subset(taxonomy,taxonomy$Kingdom!="unclassified")
counts["7751662AGE_a_G4"] <- unlist(apply(counts[c("7344072AGE1_a_G4","6897722AGE_a_G4","662586SPC_a_G4","7344072AGE1_a_G4","9457743AGE1_a_G4")][,,],1,mean))
ord = match(colnames(counts),rownames(metadata))
metadata = metadata[ord,]
#length(colnames(counts));length(union(colnames(counts),rownames(metadata)))
if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
  message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
  quit()
}
data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
counts.nl <- MRcounts(data,norm=T,log=T)
message("Filtering...")
techrep <- metadata$CC[duplicated(metadata$CC)]
xyi=c()
for (i in 1:length(techrep)){
  y <- counts.nl[,metadata[metadata$CC==techrep[i],]$ID[1]]
  x <- counts.nl[,metadata[metadata$CC==techrep[i],]$ID[2]]
  ratio <- (x+1)/(y+1)
  df <- data.frame(x=x,y=y,ratio=ratio)
  #plot(density(df$ratio,na.rm=T))
  corrxy <- cor.test(df$x,df$y)
  if(corrxy$estimate>0.77){
    minval <- mean(df$ratio)-sd(df$ratio)
    maxval <- mean(df$ratio)+sd(df$ratio)
    df$highlight <- unlist(lapply(df$ratio,function(x) if(x<minval|x>maxval){return("highlight")}else{return("normal")}))
    dz <- subset(df,df$highlight=="highlight")
    z <- abs(dz$x-dz$y)
    #xyintercept <- quantile(z,probs=0.98)
    xyintercept <- max(z)
    xyi <- c(xyi,xyintercept)
    print(c(i,corrxy$estimate,xyintercept))
    mycolours <- c("highlight" = "red", "normal" = "grey")
    ggplot(df,aes(x=x,y=y))+geom_point(aes(alpha=1/20))+
      geom_hline(yintercept=xyintercept,co)+geom_vline(xintercept=xyintercept)+
      ggtitle(i)+theme(legend.position="none",plot.title=element_text(lineheight=.8, face="bold"))
    ggsave(filename=paste(opt$out,i,'_replicates.pdf',sep=''),device="pdf")
  }
} 
rep <- metadata[metadata$CC%in%techrep,]
samplesToKeep <- rownames(metadata)[!rownames(metadata)%in%rownames(rep[rep$rep=="a",])]
data.f1 <- data[,samplesToKeep]
counts.f1.n <- MRcounts(data.f1,norm=T)
if(as.numeric(opt$shared)==0){shared <- 1}else{shared <- round(as.numeric(opt$shared)*NROW(pData(data)))}
featuresToKeep <- c()
for (i in 1:NROW(counts.f1.n)){
  x <- counts.f1.n[i,]
  x <- x[!x %in% c(0)]
  if(length(x)>=shared & mean(x)>=dep) featuresToKeep <- c(featuresToKeep,rownames(counts.f1.n)[i])
}
data.f2 <- data.f1[featuresToKeep,]
data.f <- filterData(data.f2,present=shared,depth=opt$depth)
retained.info <- sum(fData(data.f)$Size)/sum(fData(data)$Size)
fData(data.f) <- droplevels(fData(data.f))
phylum <- length(levels(as.factor(fData(data.f)$Phylum))); class <- length(levels(as.factor(fData(data.f)$Class)))
order <- length(levels(as.factor(fData(data.f)$Order))); family <- length(levels(as.factor(fData(data.f)$Family)))
genus <- length(levels(as.factor(fData(data.f)$Genus)));otu <- length(levels(as.factor(fData(data.f)$OTU)))
## Normalizing
message("Normalizing...")
p.f <- cumNormStatFast(data.f) # Calculates the percentile for which to sum counts up to and scale by.
data.f <- cumNorm(data.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(data.f)
# saving files
message("Exporting files...")
saveRDS(data.f,file=paste(opt$out,'dataF.rds',sep=''))
saveRDS(MRcounts(data.f),file=paste(opt$out,'otu_counts.rds',sep=''))
saveRDS(pData(data.f),file=paste(opt$out,'phenotype.rds',sep=''))
saveRDS(fData(data.f),file=paste(opt$out,'taxonomy.rds',sep=''))
output <- paste(paste(" - An OTU was removed when its average count was less than ",dep," and when it was present in less than the ",100*as.numeric(opt$shared),"% of the samples\n",
                      " - From ",dim(MRcounts(data))[1]," OTUs, ",dim(MRcounts(data.f))[1]," OTUs remained\n",
                      " - Information (Reads) retained: ",round(100*retained.info,2),"%; lost: ",round((100-100*retained.info),2),"%",sep=""),"\n",
                "- Taxonomic information:\n","   Phylum: ",phylum,"\tClass: ",class,"\tOrder: ",order,"\tFamily: ",family,"\tGenus: ",genus,"\tOTU: ",otu,"\n",
                "- dataF.rds -> Counts filtered and normalized - Cumulative-sum scaling normalization (Paulson et. al 2013)\n",
                "   \n** Data was successfully prepared! **\n")
message(output)
cat(output,file=paste(opt$out,'Report_data_prep.txt',sep=''), append=F)


