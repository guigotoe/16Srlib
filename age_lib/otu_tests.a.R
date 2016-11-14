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
toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","reshape2","optparse","pheatmap"))

## Options ##
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-i","--data"),type="character",default=paste(p,'age/dataF.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),type="character",default=paste(p,'age/',sep=''),#
              help="Path to output directory [default %default]"),
  make_option(c("-c","--conf"),type="character",default='Gender',#'Arcilla',#
              help="Confounder variables - separated by comma"),
  make_option(c("-v","--variable"),type="character",default='LLI',#'Salinity',#
              help="Variable of association"),
  make_option(c("-t","--shared"),type="double",default=0.15,
              help="Sample's OTU-shared percentage. 0-1; default: %default"),
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
opt$conf <- unlist(strsplit(opt$conf,','))

###### end ######
df.r <- readRDS(opt$data)
pData(df.r)$LLI <-rep('normal',NROW(pData(df.r)))
pData(df.r)$LLI[pData(df.r)$Age>90] <- "lli"
pData(df.r)$LLI <- as.factor(pData(df.r)$LLI)
#df.rr <- readRDS(opt$data)
#samplesToKeep <- which(pData(df.rr)$Gender == "female")
#df.r = df.rr[,samplesToKeep]

if (opt$level != "otu"){df <- aggTax(df.r,lvl=opt$level,norm=T)} else {df <- df.r}
### differential abundance testing -- based on OTUs-shared by 90% of the samples #### 
if(as.numeric(opt$shared)==0){th.s=0
}else{th.s <- round(as.numeric(opt$shared)*NROW(pData(df)))}

if (opt$level != "otu"){
  ## Differential abbundance test - including counfounders ###
  df.f <- filterData(df,present=th.s) # Filtration process according otu presence
  p.f <- cumNormStat(df.f,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
  cfs <- unlist(sapply(opt$conf,function(x) return(paste('pData(df.f)$',x,sep='')))) # Confounders
  
  ## FitZig Method
  
  normFactor <- normFactors(df.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
  normFactor <- log2(normFactor/median(normFactor) + 1)
  normFactor[is.na(normFactor)] <- 1
  mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df.f)[["',opt$variable,'"]])+',sep=''),
                                       paste(cfs,collapse='+'),'+normFactor')))
  settings <- zigControl(maxit=10,verbose=T)
  
  fit <- fitZig(obj=df.f,mod=mod,useCSSoffset=F,control=settings)
  #summary(calculateEffectiveSamples(fit))
  zigFit <- fit$fit
  finalMod <- fit$fit$design
  x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('factor(pData(df.f)[["',opt$variable,'"]])',sep=''),opt$variable,x,fixed=T))))
  colnames(finalMod) <- x
  x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('pData(df.f)$',sep=''),'',x,fixed=T))))
  colnames(finalMod) <- x
  colnames(zigFit$coefficients) <- colnames(finalMod)
  colnames(zigFit$stdev.unscaled) <- colnames(finalMod)
  colnames(zigFit$cov.coefficients) <- colnames(finalMod)
  colnames(zigFit$design) <- colnames(finalMod)
  k <- levels(factor(pData(df.f)[[opt$variable]]))  # states of the variable in study
  k <- unlist(sapply(k,function(x) return(paste(opt$variable,x,sep=''))))
  mk <- combn(k,2)
  contrast_list <- c()
  for (i in 1:dim(mk)[2]){contrast_list <- c(contrast_list,paste(mk[1,i],mk[2,i],sep='-'))}
  contrast.matrix <- makeContrasts(contrasts=contrast_list,levels=finalMod)
  fit2 <- contrasts.fit(zigFit,contrast.matrix)
  fit2=eBayes(fit2)
  results <- decideTests(fit2,method="separate",adjust.method="fdr",p.value=0.05) ## 1/0/-1 result DEtable 
  
  DElist <- c() 
  for (i in 1:length(colnames(fit2$coef))){
    lm <- rownames(subset(results,results[,i]!=0))
    if (length(lm)>0){
      message (colnames(fit2$coef)[i]," count talbes generation...")
      #adjustedPvalues=p.adjust(fit2$p.value[,i],method="BH")
      #foldChange=abs(fit2$coef[,i])
      #sigList=which(adjustedPvalues<=0.05)
      #sigList=sigList[order(foldChange[sigList])]
      texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=0.05,adjust="BH")
      #texp$FC <-  foldChange[rownames(texp)]
      texp$DE <- results[rownames(results)%in%rownames(texp),i]
      texp <- cbind(texp,fData(df.f)[rownames(texp),-c(1,2)])
      write.table(texp,file=paste(opt$out,opt$level,'_',colnames(fit2$coef)[i],'_',paste(opt$conf,collapse='_'),'_fitZig.txt',sep=''),
                  quote=F,sep="\t",na="NA",row.names=T)
      DElist <- c(DElist,rownames(texp))
    }
  }
  DElist.1 <- unique(DElist)
  f <- which(rownames(MRcounts(df.f))%in%DElist.1)
  dfc <- df.f[f,1:length(sampleNames(df.f))]
  
  heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(pData(dfc)$group))] 
  heatmapCols = colorRampPalette(c("#D9D9D9","#3794BF","#238B45","#D53E4F"))(4)#colorRampPalette(brewer.pal(4, "RdBu"))(50)
  for (i in levels(pData(dfc)$Gender)){
    i <- "female"
    samplesToKeep <- which(pData(dfc)$Gender==i)
    dfc.g <- dfc[,samplesToKeep]
    matc.g <- MRcounts(dfc.g,norm=T,log=T)
    pData(dfc.g) <- pData(dfc.g)[order(pData(dfc.g)$Age),]
    matc.g <- matc.g[,rownames(pData(dfc.g))]
    fontsize_col = 10 - nrow(matc.g) / 15
    fontsize_row = 10 - nrow(matc.g) / 15
    pheatmap(matc.g,col = heatmapCols,annotation_col=pData(dfc.g)[,("group"),drop=F],
             cluster_cols=F,fontsize_col=fontsize_col,fontsize_row=fontsize_row,labels_col=pData(dfc.g)[,("Age")])
    pheatmap(matc.g,col=heatmapCols,annotation_col=pData(dfc.g)[,("group"),drop=F],
             cluster_cols=T,fontsize_col=fontsize_col,fontsize_row=fontsize_row,labels_col=pData(dfc.g)[,("Age")])
  } 
  
  
  #df.f <- filterData(df,present=35) # Filtration process according otu presence
  pData(df.f) <- pData(df.f)[order(pData(df.f)$Age),]
  
  matraw <- returnAppropriateObj(df.f,norm=T,log=T)[,rownames(pData(df.f))]
  pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
           cluster_cols=F,fontsize_col=6,labels_col=pData(df.f)[,("Age")])
  pheatmap(matraw,col = heatmapCols,annotation_col=pData(df.f)[,("group"),drop=F],
           cluster_cols=T,fontsize_col=6,labels_col=pData(df.f)[,("Age")])
  fontsize_col = 10 - nrow(mat) / 15
  fontsize_row = 10 - nrow(mat) / 15

  pData(dfc) <- pData(dfc)[order(pData(dfc)$Age),]
  
  mat = returnAppropriateObj(dfc,norm=T,log=T)[,rownames(pData(dfc))]
  
  pheatmap(mat,col = heatmapCols,annotation_col=pData(dfc)[,("group"),drop=F],
           cluster_cols=F,fontsize_col=fontsize_col,fontsize_row=fontsize_row,labels_col=pData(dfc)[,("Age")])
  pheatmap(mat,col=heatmapCols,annotation_col=pData(dfc)[,("group"),drop=F],
           cluster_cols=T,fontsize_col=fontsize_col,fontsize_row=fontsize_row,labels_col=pData(dfc)[,("Age")])

  dfc <- df.f[f,1:length(sampleNames(df.f))]
  
  
  
  pie(1:12,1:12,col=brewer.pal(11,"Spectral")[9:1])
  pie(1:12,1:12,col=brewer.pal(9,"BuGn")[9:1])
  pie(1:12,1:12,col=brewer.pal(9,"Greys")[9:1])

  pie(1:10,1:10,col=colorRampPalette(c("#9BC9DF","green"))(5))
  
  c("#9BC9DF","#3794BF","#FEE08B","#D53E4F")
  counts <- MRcounts(dfc,norm=T)
  counts.l <- MRcounts(dfc,norm=T,log=T)
  pData(dfc) <- pData(dfc)[order(pData(dfc)$Age),]
  counts <- counts[,rownames(pData(dfc)[order(pData(dfc)$Age),])]
  counts.l <- counts.l[,rownames(pData(dfc)[order(pData(dfc)$Age),])]
  
  
  ### Cluster analysis
  getClusters(counts,counts.l,
              design=pData(dfc)[c("Age",opt$variable)],method=opt$clmethod,path=opt$out,prefix=paste(opt$level,'_',opt$prefix,'_',sep=""),val=opt$clval)
  ## By age#
  dfc <- df.f[f,1:length(sampleNames(df.f))]
  for (i in levels(pData(dfc)$Gender)){
    samplesToKeep <- which(pData(dfc)$Gender==i)
    dfc.g <- dfc[,samplesToKeep]
    counts <- MRcounts(dfc.g,norm=T)
    counts.l <- MRcounts(dfc.g,norm=T,log=T)
    pData(dfc.g) <- pData(dfc.g)[order(pData(dfc.g)$Age),]
    counts <- counts[,rownames(pData(dfc.g))]
    counts.l <- counts.l[,rownames(pData(dfc.g)[order(pData(dfc.g)$Age),])]
    ### Cluster analysis
    getClusters(counts,counts.l,
                design=pData(dfc.g)[c("Age",opt$variable)],method=opt$clmethod,path=opt$out,prefix=paste(opt$level,'_',opt$prefix,'_',i,'_',sep=""),val=opt$clval)
    
  } 
    
  
  
} else {
  ## ** ##
  ### Taxonomical Level analysis
  ## Differential abbundance test - including counfounders ###
  df.f <- filterData(df,present=th.s) # Filtration process according otu presence
  p.f <- cumNormStat(df.f,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
  cfs <- unlist(sapply(opt$conf,function(x) return(paste('pData(df.f)$',x,sep='')))) # Confounders
  
  ## FitZig Method
  
  normFactor <- normFactors(df.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
  normFactor <- log2(normFactor/median(normFactor) + 1)
  mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df.f)[["',opt$variable,'"]])+',sep=''),
                                       paste(cfs,collapse='+'),'+normFactor')))
  settings <- zigControl(maxit=10,verbose=T)
  
  fit <- fitZig(obj=df.f,mod=mod,useCSSoffset=F,control=settings)
  #summary(calculateEffectiveSamples(fit))
  zigFit <- fit$fit
  finalMod <- fit$fit$design
  x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('factor(pData(df.f)[["',opt$variable,'"]])',sep=''),opt$variable,x,fixed=T))))
  colnames(finalMod) <- x
  x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('pData(df.f)$',sep=''),'',x,fixed=T))))
  colnames(finalMod) <- x
  colnames(zigFit$coefficients) <- colnames(finalMod)
  colnames(zigFit$stdev.unscaled) <- colnames(finalMod)
  colnames(zigFit$cov.coefficients) <- colnames(finalMod)
  colnames(zigFit$design) <- colnames(finalMod)
  k <- levels(factor(pData(df.f)[[opt$variable]]))  # states of the variable in study
  k <- unlist(sapply(k,function(x) return(paste(opt$variable,x,sep=''))))
  mk <- combn(k,2)
  contrast_list <- c()
  for (i in 1:dim(mk)[2]){contrast_list <- c(contrast_list,paste(mk[1,i],mk[2,i],sep='-'))}
  contrast.matrix <- makeContrasts(contrasts=contrast_list,levels=finalMod)
  fit2 <- contrasts.fit(zigFit,contrast.matrix)
  fit2=eBayes(fit2)
  results <- decideTests(fit2,method="separate",adjust.method="fdr",p.value=0.05) ## 1/0/-1 result DEtable 
  DElist <- c() 
  for (i in 1:length(colnames(fit2$coef))){
    lm <- rownames(subset(results,results[,i]!=0))
    if (length(lm)>0){
      message (colnames(fit2$coef)[i]," count talbes generation...")
      #adjustedPvalues=p.adjust(fit2$p.value[,i],method="BH")
      #foldChange=abs(fit2$coef[,i])
      #sigList=which(adjustedPvalues<=0.05)
      #sigList=sigList[order(foldChange[sigList])]
      texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=0.05,adjust="BH")
      #texp$FC <-  foldChange[rownames(texp)]
      texp$DE <- results[rownames(results)%in%rownames(texp),i]
      texp <- cbind(texp,fData(df.f)[rownames(texp),-c(1,2)])
      write.table(texp,file=paste(opt$out,opt$level,'_',colnames(fit2$coef)[i],'_',paste(opt$conf,collapse='_'),'_fitZig.txt',sep=''),
                  quote=F,sep="\t",na="NA",row.names=T)
      DElist <- c(DElist,rownames(texp))
    }
  }
  DElist.1 <- unique(DElist)
  f <- which(rownames(MRcounts(df.f))%in%DElist.1)
  dfc <- df.f[f,1:length(sampleNames(df.f))]
  
  ### Cluster analysis
  df.counts <- MRcounts(dfc,norm=T)
  getClusters(df.counts,log2(df.counts+1),
              as.data.frame(pData(dfc)[[opt$variable]]),method=opt$clmethod,path=opt$out,prefix=paste(opt$level,'_',sep=""),val=3)
}

message("\nAll done!\n")

### Log Normal permutation test

#res=fitLogNormal(obj=df.f,mod=mod,useCSSoffset=F,B=10,coef=2)
#adjustedPvalues=p.adjust(res$p,method="fdr")
#foldChange=abs(res$fit$coef[,2])
#sigList=which(adjustedPvalues<=0.05)
#sigList=sigList[order(foldChange[sigList])]
#head(fData(df.f)[sigList,])

### Pair-wise comparisson using fitFeatureModel  
#comp <- combn(levels(factor(pd[[vs]])),2)
#for (i in 1:dim(comp)[2]){
#  message(comp[,i])
#  comp[,1]
#  df.f <- filterData(df,present=th.s)
#  df.f <- df.f[,which(factor(pData(df.f)[[vs]])%in%comp[,1])]
#  p.f <- cumNormStatFast(df.f) # Calculates the percentile for which to sum counts up to and scale by.
#  df.f <- cumNorm(df.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
#  s <- normFactors(df.f)
#  pd <- pData(df.f)
#  pd[[vs]] <- as.factor(pd[[vs]])
#  mod <- model.matrix(as.formula(paste("~ 1+",paste(vs,collapse="+"),sep='')),data=pd)
#  df1 <- fitFeatureModel(df.f,mod)
  #head(MRcoefs(df1))
  
#}


