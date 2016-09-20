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
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'
packages(c("metagenomeSeq","reshape2"))

###### end ######

#* input *

f <- commandArgs()[6] # paste(p,'results/dataF.rds',sep='') #
vs <- commandArgs()[7]# 'Salinity,Textura' #
#vs <- unlist(strsplit(vs,','))
cf <- commandArgs()[8]#''#,CT,NT,Ca,K,Mg,Na,CICE,Cu,S,P,Fe,Mn,Zn,B,Arcilla,Limo,Arena'
cf <- unlist(strsplit(cf,','))
th <- 0.90 
o <- commandArgs()[9] # paste(p,'results/',sep='') #
## ##
df <- readRDS(f)

### differential abundance testing -- based on OTUs shaded by 90% of the samples #### 
if(as.numeric(th)==0){th.s=0
}else{th.s <- round(as.numeric(th)*NROW(pData(df)))}

## Differential abbundance test - including counfounders ###
df.f <- filterData(df,present=th.s)
p.f <- cumNormStat(df.f,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
cfs <- unlist(sapply(cf,function(x) return(paste('pData(df.f)$',x,sep=''))))
normFactor <- normFactors(df.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
normFactor <- log2(normFactor/median(normFactor) + 1)
mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df.f)[["',vs,'"]])+',sep=''),
                                     paste(cfs,collapse='+'),'+normFactor')))
settings <- zigControl(maxit=10,verbose=T)

fit <- fitZig(obj=df.f,mod=mod,useCSSoffset=F,control=settings)
zigFit <- fit$fit
finalMod <- fit$fit$design
x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('factor(pData(df.f)[["',vs,'"]])',sep=''),vs,x,fixed=T))))
colnames(finalMod) <- x
x <- unlist(sapply(colnames(finalMod),function(x) return(gsub(paste('pData(df.f)$',sep=''),'',x,fixed=T))))
colnames(finalMod) <- x
colnames(zigFit$coefficients) <- colnames(finalMod)
colnames(zigFit$stdev.unscaled) <- colnames(finalMod)
colnames(zigFit$cov.coefficients) <- colnames(finalMod)
colnames(zigFit$design) <- colnames(finalMod)
k <- levels(factor(pData(df.f)[[vs]]))
k <- unlist(sapply(k,function(x) return(paste(vs,x,sep=''))))
c <- c(paste(combn(k,2)[1,1],combn(k,2)[2,1],sep='-'),
       paste(combn(k,2)[1,2],combn(k,2)[2,2],sep='-'),
       paste(combn(k,2)[1,3],combn(k,2)[2,3],sep='-'))
contrast.matrix <- makeContrasts(contrasts=c,levels=finalMod)
fit2 <- contrasts.fit(zigFit,contrast.matrix)
fit2=eBayes(fit2)
results <- decideTests(fit2,method="separate",adjust.method="BH",p.value=0.05)

DElist <- c() 
for (i in 1:length(colnames(fit2$coef))){
  lm <- rownames(subset(results,results[,i]!=0))
  if (length(lm)>0){
    message (colnames(fit2$coef)[i])
    adjustedPvalues=p.adjust(fit2$p.value[,i],method="fdr")
    foldChange=abs(fit2$coef[,i])
    #sigList=which(adjustedPvalues<=0.05)
    #sigList=sigList[order(foldChange[sigList])]
    texp <- topTable(fit2,coef=i,number=10000,p.value=0.05,adjust="BH")
    texp$FC <-  foldChange[rownames(texp)]
    texp <- cbind(texp,fData(df.f)[rownames(texp),-c(1,2)])
    write.table(texp,file=paste(o,colnames(fit2$coef)[i],'_',paste(cf,collapse='_'),'_fitZig.txt',sep=''),
                quote=F,sep="\t",na="NA",row.names=T)
    DElist <- c(DElist,rownames(texp))
  }
}
DElist.1 <- unique(DElist)
f <- which(rownames(MRcounts(df.f))%in%DElist.1)
dfc <- df.f[f,1:length(sampleNames(df.f))]

### Cluster analysis

getClusters(MRcounts(dfc,norm=T),log2(MRcounts(dfc,norm=T)+1),
            as.data.frame(pData(dfc)[[vs]]),method="K",path=o,prefix="",val=4)

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


