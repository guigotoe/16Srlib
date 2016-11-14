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
# Rscript otu_test.R /path/dataF.rds association_variable co-variants/path/outfolder/
# Rscript otu_test.R ~/16Srlib_test/results/dataF.rds Salinity Limo,Arena ~/16Srlib_test/results/
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

packages(c("metagenomeSeq","reshape2","optparse","pheatmap","vegan","clusterSim","Rlof"))

## Options ##
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-i","--data"),type="character",default=paste(p,'age/dataF.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),type="character",default=paste(p,'age/',sep=''),#
              help="Path to output directory [default %default]"),
  make_option(c("-t","--shared"),type="double",default=100,
              help="Sample's OTU-shared 0-n; default: %default"),
  make_option(c("-l","--level"),type="character",default="otu",
              help="Taxonomical level of the analysis (otu,Genus,Family,Order,Class,Phylum). default: %default"),
  make_option(c("-m","--clmethod"),type="character",default='P',
              help="Clusterization method.(PAM,P,K,Km)\nUse K-means clustering to define K feature sets \n
              ; default: %default"),
  make_option(c("-n","--clval"),type="double",default=80,
              help="Number of clusters (K method) or \
              cut the hierarchically clustered tree at -n percent height of the tree (P method).\n
              default: %default"),
  make_option(c("-p","--prefix"),type="character",default='LLI',#'Arcilla',#
              help="files prefix")
  )
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not data file specified"))}
###### end ######
###### reading files ######
df.r <- readRDS(opt$data)
pData(df.r) <- pData(df.r)[order(pData(df.r)$Age),]  # soring phenotype according Age
df.r <- df.r[,rownames(pData(df.r))]                 # soring samples according Age
pData(df.r)$LLI <-rep('normal',NROW(pData(df.r)))
pData(df.r)$LLI[pData(df.r)$Age>90] <- "lli"
pData(df.r)$LLI <- as.factor(pData(df.r)$LLI)
pData(df.r)$group <- as.factor(unlist(apply(pData(df.r),1,function(x){
  if (as.numeric(x[["Age"]]) < 40) {return("G1")
  }else if (as.numeric(x[["Age"]]) > 39 & as.numeric(x[["Age"]]) < 60){ return("G2")
  }else if (as.numeric(x[["Age"]]) > 59 & as.numeric(x[["Age"]]) < 80){ return("G3")
  }else if (as.numeric(x[["Age"]]) > 79) return("G4")
})))
###### filtrating ######
df.f <- filterData(df.r,present=opt$shared) # Filtration process according otu presence

## heatmap of the OTU presence ##
heatmapColColors <- brewer.pal(9, "Set1")[as.integer(factor(pData(df.f)$group))] 
heatmapCols <- colorRampPalette(c("#D9D9D9","#3794BF","#238B45","#D53E4F"))(4)
fontsize <- 6
matraw <- returnAppropriateObj(df.f,norm=T,log=T)
pdf(file=paste(opt$out,"heatmap_all_clF.pdf",sep=""),width=9,height=7,onefile=FALSE)
pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
         cluster_cols=F,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.f)[,("Age")])
dev.off()
pdf(file=paste(opt$out,"heatmap_all_clT.pdf",sep=""),width=9,height=7,onefile=FALSE)
pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
         cluster_cols=T,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.f)[,("Age")])
dev.off()
for (i in levels(pData(df.f)$Gender)){
  samplesToKeep <- which(pData(df.f)$Gender==i)
  dfc.g <- df.f[,samplesToKeep]
  matc.g <- MRcounts(dfc.g,norm=T,log=T)
  pdf(file=paste(opt$out,"heatmap_",i,"_clF.pdf",sep=""),width=9,height=7,onefile=FALSE)
  pheatmap(matc.g,col = heatmapCols,annotation_col=pData(dfc.g)[,("group"),drop=F],
           cluster_cols=F,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(dfc.g)[,("Age")])
  dev.off()
  pdf(file=paste(opt$out,"heatmap_",i,"_clT.pdf",sep=""),width=9,height=7,onefile=FALSE)
  pheatmap(matc.g,col=heatmapCols,annotation_col=pData(dfc.g)[,("group"),drop=F],
           cluster_cols=T,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(dfc.g)[,("Age")])
  dev.off()
} 

### geting males and females data ##
samplesToKeep <- which(pData(df.f)$Gender=="male")
males <- df.f[,samplesToKeep]
samplesToKeep <- which(pData(df.f)$Gender=="female")
females <- df.f[,samplesToKeep]
colorRampPalette( brewer.pal(9, "Set1"))(6)

Nets <- list()
phenotype.f <- as.data.frame(matrix(data=NA,nrow=0,ncol=length(colnames(pData(females)))))
colnames(phenotype.f) <- colnames(pData(females))
j="males"
i <- "G4"
df <- males
### by age group ###
for (i in pData(df)$group){
  samplesToKeepByGroup <- which(pData(df)$group==i)
  df.g <- df[,samplesToKeepByGroup]
  ## remove outliers using Local Outlier Factor (LOF) approach##
  dfc <- t(MRcounts(df.g,norm=T,log=T))
  distmat <- vegdist(dfc,"bray")
  nmds <- metaMDS(distmat,trace=F)
  ordiplot(nmds,type="p",display="sites")
  sc <- scores(nmds)
  outliers <- which(lof(sc,nrow(dfc)/3)>2)#
  plot(lof(sc,nrow(dfc)/3))
  if(length(outliers)!=0) df.g <- df.g[,-outliers]
  phenotype.f <- rbind(phenotype.f,pData(df.g))
  ## end outliers##
  
  matf.g <- MRcounts(df.g,norm=T,log=T)
  #pdf(file=paste(opt$out,"heatmap_",i,"_clF.pdf",sep=""),width=9,height=7,onefile=FALSE)
  pheatmap(matf.g,col = heatmapCols,annotation_col=pData(df.g)[,("group"),drop=F],
           cluster_cols=F,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.g)[,("Age")])
  #dev.off()
  #pdf(file=paste(opt$out,"heatmap_",i,"_clT.pdf",sep=""),width=9,height=7,onefile=FALSE)
  pheatmap(matf.g,col=heatmapCols,annotation_col=pData(df.g)[,("group"),drop=F],
           cluster_cols=T,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.g)[,("Age")])
  #dev.off()
  
  ## getting clusters##
  data <- matf.g
  #sample_dist <- vegdist(t(data),"bray")
  data.dist=dist.JSD(data)
  #data.cluster=pam.clustering(data.dist, k=3)
  #idx <- index.G1(data,data.cluster,d=data.dist,centrotypes="medoids")
  nclusters=NULL
  for (k in 1:round(ncol(matf.g)/2)) { 
    if (k==1) {
      nclusters[k]=NA 
    } else {
      data.cluster_temp=pam.clustering(data.dist, k=k)
      nclusters[k]=index.G1(t(data),data.cluster_temp,d=data.dist,centrotypes="medoids")
    }
  }
  #nclusters
  #pdf(paste(opt$out,"OptClusterNum_CHindex.pdf",sep=''),width=8, height=5)
  plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
  #dev.off()
  #val <- which(nclusters==max(nclusters[-1]))[1]
  nc <- which(nclusters==max(nclusters[-1]))
  data.cluster=pam.clustering(data.dist, k=nc)
  s <- data.frame(id=colnames(data),clusters=data.cluster)
  s <- s[with(s, order(clusters)),]
  df.g <- df.g[,as.character(s$id)]
  matf.gc <- MRcounts(df.g,norm=T,log=T)
  pheatmap(matf.gc,main=paste(i,j,sep=" "),co=heatmapCols,annotation_col=pData(df.g)[,("group"),drop=F],
           cluster_cols=F,fontsize_col=fontsize,fontsize_row=fontsize,labels_col=pData(df.g)[,("Age")])
  ## coocurrence for each cluster ##
  net <- net.mb(df.g)
  title <- paste(j,i,sep="_")
  Nets[[title]] <- net
  pdf(paste(opt$out,title,"_net.pdf",sep=''),width=8, height=5)
  plot(net$graph, layout=net$layout, vertex.size=net$vsize, vertex.label=NA, main=paste(j,i,sep=" "),edge.width=2)#abs(E(net$graph)$weight)*40
  dev.off()  
  k <- 1
  df.gc <- df.g[,as.character(s[s$clusters==k,]$id)]
  net <- net.mb(df.gc)
  title <- paste(j,i,"cluster",k,sep="_")
  Nets[[title]] <- net
  pdf(paste(opt$out,title,"_net.pdf",sep=''),width=8, height=5)
  plot(net$graph, layout=net$layout, vertex.size=net$vsize, vertex.label=NA, main=paste(j,i,sep=" "),edge.width=2)
  dev.off() 
  }
}

