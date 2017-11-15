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
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/age_lib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse","ggplot2","outliers"))

## Options ##
pa = '/home/torres/ikmb_storage/Metagenome/16s/2017/4_2017/'#
ra <- '/home/torres/Documents/Projects/Metagenome/2017_results/9_2017/'
#pm <- '/home/torres/ikmb_storage/Mangrove/16Sfa/2017/mothur/silva_result/'
#rm <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2017/'

option_list <- list(
  make_option(c("-c","--counts"),action="store",type="character",default=paste(pa,'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.subsample.opti_mcc.unique_list.shared',sep=''),#'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared'
              help="Path to input counts file"),
  make_option(c("-x","--taxonomy"),action="store",type="character",default=paste(pa,'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.subsample.opti_mcc.unique_list.0.03.cons.taxonomy',sep=''),# '16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy'
              help="Path to input taxonomy file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(pa,'metadata.txt',sep=''),
              help="Path to input metadata file"),
  make_option(c("-f","--fasta"),action="store",type="character",default=paste(pa,'OTUs.fa',sep=''),
              help="Path to OTUs fasta file"),
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
#counts <- subset(counts,taxonomy$Kingdom!="unclassified"& taxonomy$Phylum!="unclassified")
#taxonomy <- subset(taxonomy,taxonomy$Kingdom!="unclassified"& taxonomy$Phylum!="unclassified")
#counts["7751662AGE_a_G4"] <- unlist(apply(counts[c("7344072AGE1_a_G4","6897722AGE_a_G4","662586SPC_a_G4","7344072AGE1_a_G4","9457743AGE1_a_G4")][,,],1,mean))
ord = match(colnames(counts),rownames(metadata))
metadata = metadata[ord,]
length(colnames(counts));length(union(colnames(counts),rownames(metadata)))
if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
  message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
  quit()
}
dim(counts)
dim(taxonomy)
data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))

message("Filtering...")

## technical replicates ?? ##

if(opt$replicates){
  techrep <- replicas.analysis(data,opt$out)
  data.f1 <- techrep$data
  otu.coverage <- techrep$coverage
  otu.coverage <- round(2^(mean(techrep$interceptcoords)))
  counts.f1.n <- MRcounts(data.f1,norm=T)
  #counts.f1.n <- counts.nl[,which(samplesToKeep%in%colnames(counts.nl))]
}else {
  otu.coverage <- opt$depth
  data.f1 <- data
  counts.f1.n <- MRcounts(data.f1,norm=T)
}
### # cutoff at %shared% counts in at least X libraries.
if (is.integer(opt$libraries) & opt$libraries>1){
  if(as.numeric(opt$shared)==0){shared <- 1}else{shared <- round(as.numeric(opt$shared)*NROW(pData(data)))}
  ## Identifying which OTUs satisfying only coverage threshold
  featuresToKeep.a <- c()
  for (i in 1:NROW(counts.f1.n)){
    x <- counts.f1.n[i,]
    x <- x[!x %in% c(0)]
    z <- x>otu.coverage
    if(length(x)>=1 & length(z[z==TRUE]) >=1) featuresToKeep.a <- c(featuresToKeep.a,rownames(counts.f1.n)[i])
  }
  ## Identifying which OTUs satisfying presence and coverage thresholds
  featuresToKeep.b <- c()
  for (i in 1:NROW(counts.f1.n)){
    x <- counts.f1.n[i,]
    x <- x[!x %in% c(0)]
    z <- x>otu.coverage
    if(length(x)>=shared & length(z[z==TRUE]) >= opt$libraries) featuresToKeep.b <- c(featuresToKeep.b,rownames(counts.f1.n)[i])
  }
  D2 <- TRUE
}else{
  ## adjusting shared (samples sharing otus)
  #opt$shared <- 0.2
  if(as.numeric(opt$shared)==0){shared <- 1}else{shared <- round(as.numeric(opt$shared)*NROW(pData(data)))}
  
  ## Identifying which OTUs satisfying only coverage threshold
  featuresToKeep.a <- c()
  for (i in 1:NROW(counts.f1.n)){
    x <- counts.f1.n[i,]
    x <- x[!x %in% c(0)]
    if(length(x)>=1 & mean(x)>=otu.coverage) featuresToKeep.a <- c(featuresToKeep.a,rownames(counts.f1.n)[i])
  }
  ## Identifying which OTUs satisfying presence and coverage thresholds
  featuresToKeep.b <- c()
  for (i in 1:NROW(counts.f1.n)){
    x <- counts.f1.n[i,]
    x <- x[!x %in% c(0)]
    if(length(x)>=shared & mean(x)>=otu.coverage) featuresToKeep.b <- c(featuresToKeep.b,rownames(counts.f1.n)[i])
  }
  D2 <- FALSE
}
data.fa <- data.f1[featuresToKeep.a,]
data.a <- filterData(data.fa,present=1,depth=opt$depth)
data.fb <- data.f1[featuresToKeep.b,]
data.b <- filterData(data.fb,present=shared,depth=opt$depth)
packages("Biostrings")
OTUs.fa <- readDNAStringSet(opt$fasta)
for (i in 1:2){
  if(i==1){
    if(D2==TRUE) {data.f <- data.a;fname <- "c_l";shared.prc <- paste("in at least 1 libraries\n",sep='');avg <- ''}
    if(D2==FALSE) {data.f <- data.a;fname <- "c";shared.prc <- paste("and when it was present in at least 1 sample\n",sep='') ;avg <- 'average '}
  }
  if (i==2){
    if(D2==TRUE) {data.f <- data.b;fname <- paste("cp_l_",opt$shared,sep='');shared.prc <- paste("in at least ",opt$libraries," libraries\n",sep='');avg <- ''}
    if(D2==FALSE) {data.f <- data.b;fname <- paste("cp_",opt$shared,sep='');shared.prc <- paste("more than ",shared,' (',100*as.numeric(opt$shared),'%)'," samples\n",sep='');avg <- 'average '}
  }
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
  # Fasta files
  if(names(OTUs.fa)[1]!=as.character(fData(data.f)$OTU)[1]) names(OTUs.fa) <- gsub("Otu","Otu0",names(OTUs.fa),fixed=T)
  selected.fa <- OTUs.fa[names(OTUs.fa)%in%as.character(fData(data.f)$OTU)]
  writeXStringSet(selected.fa,paste(opt$out,'OTU_',fname,'.fa',sep=''),format="fasta")
  ###
  saveRDS(data.f,file=paste(opt$out,'dataF',fname,'.rds',sep=''))
  saveRDS(MRcounts(data.f,norm=T),file=paste(opt$out,fname,'_','otu_counts.norm.rds',sep=''))
  saveRDS(pData(data.f),file=paste(opt$out,fname,'_','phenotype.rds',sep=''))
  saveRDS(fData(data.f),file=paste(opt$out,fname,'_','taxonomy.rds',sep=''))
  output <- paste(paste("Filtred mode ",fname,":\n An OTU was kept when its ",avg,"count was grater than ",otu.coverage," reads ",shared.prc,
                        " - From ",dim(MRcounts(data))[1]," OTUs, ",dim(MRcounts(data.f))[1]," OTUs remained\n",
                        " - Information: ",sum(fData(data)[2])," Total reads;"," retained: ",round(100*retained.info,2),"%; lost: ",round((100-100*retained.info),2),"%",sep=""),"\n",
                  "- Taxonomic information:\n","   Phylum: ",phylum,"\tClass: ",class,"\tOrder: ",order,"\tFamily: ",family,"\tGenus: ",genus,"\tOTU: ",otu,"\n",
                  "- dataF.rds -> Counts filtered and normalized - Cumulative-sum scaling normalization (Paulson et. al 2013)\n",
                  "   \n** Data was successfully prepared! **\n")
  message(output)
  cat(output,file=paste(opt$out,'Report_data_prep.txt',sep=''),append=T)
  
  
  ## ##
  
}
remove(counts.f1.n)

##### tech replicates ###
#techrep <- metadata$CC[duplicated(metadata$CC)]
#packages(c('gridExtra'))
#xyi=c()
#plots <- list()
#for (i in seq(length(techrep))){
  #vec <- 1:nrow(metadata[metadata$CC==techrep[i],])
  #combs <- combn(vec,2)
  #for (j in 1:ncol(combs)){
  #  df.s <- metadata[metadata$CC==techrep[i],][combs[,j],]
  #  y <- counts.nl[,df.s$ID[1]]
  #  x <- counts.nl[,df.s$ID[2]]
  #  df <- data.frame(x=x,y=y)
  #  df$ratio <- unlist(apply(df,1,function(x){
  #    if(is.finite(x[1]/x[2])){if(x[1]/x[2]==0)return(x[2]) else return(x[1]/x[2])
  #    }else if(is.nan(x[1]/x[2])){return(1)
  #    }else if(is.infinite(x[1]/x[2])) return(x[1])
  #  }))
  #  #plot(density(df$ratio))
  #  corrxy <- cor.test(df$x,df$y)
  #  if(corrxy$estimate>0.8){
  #    minval <- mean(df$ratio)-sd(df$ratio)
  #    maxval <- mean(df$ratio)+sd(df$ratio)
  #    df$highlight <- unlist(lapply(df$ratio,function(x) if(x!=1){return("highlight")}else{return("normal")}))
  #    dz <- subset(df,df$highlight=="highlight")
  #    z <- abs(dz$x-dz$y) # how far is one value to another.
  #    xyintercept <- quantile(z,probs=0.95)
  #    #hist(z)
  #    #abline(v=xyintercept)
  #    #xyintercept <- max(z)
  #    xyi <- c(xyi,xyintercept)
  #    df$keept <- unlist(apply(df,1,function(x) if(as.numeric(x[1])>xyintercept&as.numeric(x[2])>xyintercept){return("keept")}else{return("drop")}))
  #    print(c(paste(i,j,sep="_"),corrxy$estimate,xyintercept,'OTUs_retained'=length(df$keept[df$keept=='keept'])))
  #    mycolours <- c("keept" = "red", "drop" = "black")
  #    dz2 <- subset(df,df$highlight=="highlight")
  #    fig <- ggplot(dz2,aes(x=x,y=y))+geom_point(aes(alpha=1/100,colour=keept))+
  #      geom_hline(yintercept=xyintercept,co)+geom_vline(xintercept=xyintercept)+scale_colour_manual(values=mycolours)+
  #      ggtitle(paste(i,'_',j," Pcor: ",round(corrxy$estimate,2),' Int.val(log2): ',round(xyintercept,2),' Int.val(counts): ',round(2^xyintercept,2),sep=""))+
  #      theme(legend.position="none",plot.title=element_text(lineheight=.8, face="bold"))
  #    plots[[paste(i,j,sep="_")]] <- fig
  #    #ggsave(filename=paste(opt$out,i,'_replicates.pdf',sep=''),device="pdf")
  #  }
  #}
 # pdf(paste(opt$out,'TechRep_scatter.pdf',sep=''),onefile=T,width=9,height=7)
 # for (q in seq(length(plots))){do.call("grid.arrange",plots[q])}
#  dev.off()
#}
#plotDensities(xyi)
##threshold <- quantile(xyi,probs=0.95)
#shapiro.test(xyi)
#hist(xyi)
#abline(v=threshold)
#abline(v=(mean(xyi)+sd(xyi)))0
#round(2^mean(threshold))
#summary(xyi)
#round(2^threshold)
#round(2^mean(xyi))
#round(2^(mean(xyi)+sd(xyi)))
#otu.coverage <- round(2^(mean(xyi)+sd(xyi))) # threshold implemented as mean+sd. Stringent enough infered from data.
###otu.coverage <- round(2^mean(5.5)) ## manual adjust
### removing replications
#rep <- metadata[metadata$CC%in%techrep,]
##samplesToKeep <- rownames(metadata)[!rownames(metadata)%in%rownames(rep[rep$rep=="a",])]
#samplesToKeep <- rownames(metadata[is.na(metadata$techrep),])
#data.f1 <- data[,samplesToKeep]
#counts.f1.n <- MRcounts(data.f1,norm=T)
#counts.f1.n <- counts.nl[,which(samplesToKeep%in%colnames(counts.nl))]
##otu.coverage <- 6 ##*** because I did prev. manually and is log=T then i should use xyi value: manualy placed as 6


