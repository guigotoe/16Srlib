####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: November 2017
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
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","reshape2","optparse","vegan"))

## Options ##
p <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-i","--data"),type="character",default=paste(p,'dataFcp_l_0.1.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),type="character",default=paste(p,'fitZig_noCovar/',sep=''),#
              help="Path to output directory [default %default]"),
  make_option(c("-c","--conf"),type="character",default='',#
              help="Confounder variables - separated by comma"),
  make_option(c("-v","--variable"),type="character",default='Salinity',#
              help="Variable of association"),
  make_option(c("-p","--pval"),type="double",default=0.05,
              help="Significance of adjpval.default: %default"),
  make_option(c("-t","--shared"),type="double",default=0.05,
              help="Sample's OTU-shared percentage. 0-1; default: %default"),
  make_option(c("-l","--level"),type="character",default="Genus",
              help="Taxonomical level of the analysis (OTU,Genus,Family,Order,Class,Phylum). default: %default"),
  make_option(c("-m","--clmethod"),type="character",default='P',
              help="Clusterization method.(PAM,P,K,Km)\nUse K-means clustering to define K feature sets \n
              ; default: %default"),
  make_option(c("-n","--clval"),type="double",default=70,
              help="Number of clusters (K method) or \
              cut the hierarchically clustered tree at -n percent height of the tree (P method).\n
              default: %default")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not data file specified"))}
opt$conf <- unlist(strsplit(opt$conf,','))
if(dir.exists(opt$out)){message('Out-folder already exist, files will be overwritten')
}else dir.create(opt$out,showWarnings=F)

refvar <- NA ## reference variable... to make the plots according to this varialbe
#refvar <- 'ID_ref'

###### end ######
df.r <- readRDS(opt$data)
pData(df.r) <- pData(df.r)[,colSums(is.na(pData(df.r)))==0] ## remove coluns with NAs

## use of refvar ##
refvar <- 'ID_ref'
pData(df.r)[[refvar]] <- factor(pData(df.r)[[refvar]],levels=c('low','med','high'))
#### end refvar ##

if (opt$level != "OTU"){df <- aggTax(df.r,lvl=opt$level)
#counts <- t(decostand(t(MRcounts(df)),method='hellinger')) ##Helinger transformation is useless.. doesnt works for this.. just with cumSum normalization
#df <- newMRexperiment(counts, phenoData=AnnotatedDataFrame(pData(df)),featureData=AnnotatedDataFrame(fData(df)))
} else {df <- df.r}
write.table(round(MRcounts(df,norm=T),2),file=paste(opt$out,opt$level,'_AllNormCounts.txt',sep=''),
            quote=F,sep="\t",na="NA",row.names=T)

### differential abundance testing -- based on OTUs-shared by 90% of the samples #### 
if(as.numeric(opt$shared)==0){th.s=0
}else{th.s <- round(as.numeric(opt$shared)*NROW(fData(df)))}
fit <- NA
if (!is.na(opt$level)){
  ## Differential abbundance test - including counfounders ###
  #df.f <- filterData(df,present=th.s) # Filtration process according otu presence
  p.f <- cumNormStat(df,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
  df <- cumNorm(df,p=p.f)
  cfs <- unlist(sapply(opt$conf,function(x) return(paste('pData(df)$',x,sep='')))) # Confounders
  
  ## FitZig Method
  normFactor <- normFactors(df)  # Calculates each column's quantile and calculates the sum up to and including p quantile
  #normFactor <- log2(normFactor/median(normFactor) + 1)
  #normFactor[is.na(normFactor)] <- 1
  mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df)[["',opt$variable,'"]])+',sep=''),
                                       paste(cfs,collapse='+'),'+normFactor')))
  settings <- zigControl(maxit=10,verbose=T)
  fit <- NULL
  try(fit <- fitZig(obj=df,mod=mod,control=settings,useMixedModel=T))
  if(length(fit)==0){
    fit <- fitZig(obj=df,mod=mod,control=settings,useMixedModel=F)
  }
  #summary(calculateEffectiveSamples(fit))
  zigFit <- fit$fit
  finalMod <- fit$fit$design
  x <- unlist(lapply(colnames(finalMod),function(x) return(gsub(paste('factor(pData(df)[["',opt$variable,'"]])',sep=''),opt$variable,x,fixed=T))))
  colnames(finalMod) <- x
  x <- unlist(lapply(colnames(finalMod),function(x) return(gsub('pData(df)$','',x,fixed=T))))
  colnames(finalMod) <- x
  colnames(zigFit$coefficients) <- colnames(finalMod)
  colnames(zigFit$stdev.unscaled) <- colnames(finalMod)
  colnames(zigFit$cov.coefficients) <- colnames(finalMod)
  colnames(zigFit$design) <- colnames(finalMod)
  k <- levels(factor(pData(df)[[opt$variable]]))  # states of the variable in study
  k <- unlist(sapply(k,function(x) return(paste(opt$variable,x,sep=''))))
  mk <- combn(k,2)
  contrast_list <- c()
  for (i in 1:dim(mk)[2]){contrast_list <- c(contrast_list,paste(mk[1,i],mk[2,i],sep='-'))}
  contrast.matrix <- makeContrasts(contrasts=contrast_list,levels=finalMod)
  fit2 <- contrasts.fit(zigFit,contrast.matrix)
  fit2=eBayes(fit2)
  results <- decideTests(fit2,method="separate",adjust.method="fdr",p.value=0.01,lfc=0.5) ## 1/0/-1 result DEtable 
  DElist <- c() 
  for (i in 1:length(colnames(fit2$coef))){
    lm <- rownames(subset(results,results[,i]!=0))
    if (length(lm)>0){
      message (colnames(fit2$coef)[i]," count talbes generation...")
      texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=opt$pval,adjust="BH")
      if (nrow(texp)>0){
        texp$DE <- results[rownames(results)%in%rownames(texp),i]
        texp <- cbind(texp,fData(df)[rownames(texp),-c(1,2)])
        write.table(texp,file=paste(opt$out,opt$level,'vs',colnames(fit2$coef)[i],'_conf_',paste(opt$conf,collapse='_'),'_fitZig.txt',sep=''),
                    quote=F,sep="\t",na="NA",row.names=T)
        DElist <- c(DElist,rownames(texp))
      }
    }else{
      message(paste(colnames(fit2$coef)[i],' has no significant elements at ',opt$pval,' threshold'))
      #texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=opt$pval,adjust="BH")
      #write.table(texp,file=paste(opt$out,opt$level,'vs',colnames(fit2$coef)[i],'_conf_',paste(opt$conf,collapse='_'),'_fitZig.txt',sep=''),
      #            quote=F,sep="\t",na="NA",row.names=T)
    }
  }
  DElist.1 <- unique(DElist)
  f <- which(rownames(MRcounts(df))%in%DElist.1)
  dfc <- df[f,1:length(sampleNames(df))]
  dfexport_counts <- setNames(data.frame(round(MRcounts(dfc,norm=T),2)),colnames(MRcounts(dfc)))
  if(opt$level=='OTU'){
    dfexport_alltax <- lapply(rownames(dfexport_counts),function(x){
      fData(df.r)[fData(df.r)['OTU']==x,][1,3:8]
    })
  }else{
    dfexport_alltax <- lapply(rownames(dfexport_counts),function(x){
      fData(df.r)[fData(df.r)[opt$level]==x,][1,3:which(colnames(fData(df.r))==opt$level)]
    })
  }
    
  dfexport_alltax.df <- do.call(rbind.data.frame,dfexport_alltax)
  dfexport <- cbind(dfexport_counts,dfexport_alltax.df)
  write.table(dfexport,file=paste(opt$out,opt$level,'DifAbund_conf_',paste(opt$conf,collapse='_'),'_CountsTaxonomy.txt',sep=''),
              quote=F,sep="\t",na="NA",row.names=T)
  ### Cluster analysis
  #source(toolbox)
  if(nrow(MRcounts(dfc))>0){
    if(!is.na(refvar)) {design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,refvar])),c('ID',opt$variable))
    }else design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,opt$variable])),c('ID',opt$variable))
    design <- design[with(design,order(design[,2])),]
    rownames(design) <- design$ID
    getClusters( MRcounts(dfc,norm=T), MRcounts(dfc,norm=T,log=T),
                 design=design,method=opt$clmethod,path=opt$out,
                 prefix=paste(opt$level,'Sign_conf_',paste(opt$conf,collapse='_'),'_',sep=""),val=opt$clval,
                 cellwidth=20,cellheight=8.3,text_size=8.7,height=13,width=7)
  }
  ## ploting all genders..
  #if(!is.na(refvar)) {design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,refvar])),c('ID',opt$variable))
  #}else design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,opt$variable])),c('ID',opt$variable))
  #design <- design[with(design,order(design[,2])),]
  #rownames(design) <- design$ID
  #getClusters(MRcounts(df,norm=T), MRcounts(df,norm=T,log=T),
  #             design=design,method=opt$clmethod,path=opt$out,prefix=paste(opt$level,'All_',sep=""),val=50)
  
} 
message("\nAll done!\n")
fitlognorm==F
if (fitlognorm){
  ## Differential abbundance test - including counfounders ###
  ## filtrating sample ##
  pData(df)
  df.x <- df[,-which(pData(df)$ID_ref=='med')]
  p.f <- cumNormStat(df.x,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
  df.x <- cumNorm(df.x,p=p.f)
  cfs <- unlist(sapply(opt$conf,function(x) return(paste(x,sep='')))) # Confounders
  
  ## FitLogNormal Method
  normFactor <- normFactors(df.x)  # Calculates each column's quantile and calculates the sum up to and including p quantile
  normFactor <- log2(normFactor/median(normFactor) + 1)
  normFactor[is.na(normFactor)] <- 1
  mod <- model.matrix(as.formula(paste("~",paste('1+',opt$variable,'+',sep=''),
                                       paste(cfs,collapse='+'),'+normFactor')),data=pData(df.x))
  
  mod <- model.matrix(as.formula(paste("~",paste('1+',opt$variable,sep=''))),data=pData(df.x))
  
  res <- fitLogNormal(obj=df.x,mod=mod)
  MRcoefs(res)
  res$fit
} 

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


