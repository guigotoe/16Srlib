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
#toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
packages(c("metagenomeSeq","reshape2","vegan","ggplot2","optparse","matrixStats"))

## Options ##
#p <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib_test/age/'
p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/age/'

option_list <- list(
  make_option(c("-i","--data"),action="store",type="character",default=paste(p,'dataF.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'beta',sep=''),#"./",#
             help="Path to output directory [default %default]"),
  make_option(c("-e","--exploratory"),action="store_true",default="NA",
              help="Perform exploratory analysis"),
  make_option(c("-m","--model"),action="store_true",default=FALSE,
              help="Build constrained model based on AIC selection criterion"),
  make_option(c("-am","--assess_model"),action="store",type="character",default=NA,#NULL,
              help="Model's terms assessed by permutation tests; for new model: vars,separated,by,comma"),
  make_option(c("-C","--constraints"),action="store",type="character",default=NA,
              help="Set of constraints used by -coa: vars,separated,by,comma"),
  make_option(c("-coa","--constrained_analysis"),action="store_true",default=FALSE,
              help="Perform constrained ordination analysis using -C constraints"),
  make_option(c("-a","--factors"),action="store",type="character",default=NULL,
              help="Set of factors used by -b: vars,separated,by,comma"),
  make_option(c("-b","--beta"),action="store_true",default=FALSE,
              help="Perform beta diversty between -a variable classes"),
  make_option(c("-f","--filter"),action="store",type="double",default=0.3,
              help="Percentile as threshold of low abundant features ")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not file specified"))}
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

#### Preparing the input data ####
data <- readRDS(opt$data)
mat <- MRcounts(data,norm=T)
otuMeans <- rowMeans(mat)
otuMeans <- otuMeans[order(otuMeans,decreasing=F)]
totrim <- names(otuMeans)[which(otuMeans<=quantile(otuMeans,probs=opt$filter))]
totrim <- match(totrim,rownames(data))
dataTrimed <- data[-totrim,]

dfc <- t(MRcounts(dataTrimed,norm=T,log=T))
#ptrim <- pData(dataTrimed)
#ptrim$group <- as.factor(unlist(apply(ptrim,1,function(x){
#  if (as.numeric(x[["Age"]]) <= 40) {return("G1")
#  }else if (as.numeric(x[["Age"]]) > 40 & as.numeric(x[["Age"]]) <= 60){ return("G2")
#  }else if (as.numeric(x[["Age"]]) > 60 & as.numeric(x[["Age"]]) <= 80){ return("G3")
#  }else if (as.numeric(x[["Age"]]) > 80) return("G4")
#})))
pData(dataTrimed)$LLI <-rep('normal',NROW(pData(dataTrimed)))
pData(dataTrimed)$LLI[pData(dataTrimed)$Age>90] <- "lli"
pData(dataTrimed)$LLI <- as.factor(pData(dataTrimed)$LLI)
q <- pData(dataTrimed)[, !sapply(pData(dataTrimed), is.factor),drop=F]
c <- pData(dataTrimed)[, sapply(pData(dataTrimed), is.factor),drop=F]
#taxa <- fData(df)[,which(colnames(fData(df))%in%c("Phylum","Class","Order","Family"))] 
#rankindex(scale(x),t(MRcounts(df,norm=T)),c("euc","man","bray","jac","kul"))

###### end ######

if (opt$exploratory==T){
  pdf(paste(opt$out,"e01_data_filtering.pdf",sep=''),width=8, height=5)
  plot(otuMeans)
  abline(h=quantile(otuMeans,probs=opt$filter))
  dev.off()
  #dist <- vegdist(dfc)
  # 1) Direction of the gradient: The arrow points of the direction of most rapid change 
  #    in the environmental variable.
  # 2) Strength of the gradient: The length of the arrow is proportional to the correlation 
  #    between ordination and environmental variable.
  
  distmat <- vegdist(dfc,"bray")
  ord = cmdscale(distmat,k = max(2))
  xl = paste("MDS component:",1)
  yl = paste("MDS component:",2)
  plot(ord[,2],ylab=yl,xlab=xl)
  
  nmds <- metaMDS(dfc,trace=F)
  pdf(paste(opt$out,"e02_NMDS.pdf",sep=''),width=8, height=5)
  ordiplot(nmds,type="p",display="sites")
  #with(c,ordiellipse(nmds,group,kind="se",conf=0.95))
  #with(c,ordispider(nmds,group,col="blue",label="T"))
  #with(c,ordihull(nmds,group,col="blue",lty=2))
  dev.off()
  ef <- envfit(nmds,q[,-which(colnames(q)=="libsize"),drop=F],permu=999)
  sink(file=paste(opt$out,"e_Envfitting2NMDS.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e03_EnvNMDS.pdf",sep=''),width=8, height=5)
  plot(nmds,type="n")
  with(c, points(nmds, display = "sites", col = cols[group],pch = 16,bg=cols[group]))
  plot(ef,p.max=0.05)
  dev.off()
  
  pca <- rda(dfc,scale=TRUE)
  pdf(paste(opt$out,"e04_PCA.pdf",sep=''),width=8, height=5)
  plot(pca)
  #text(pca,type="p",display="sites")
  dev.off()
  ef <- envfit(pca,q[,-which(colnames(q)=="libsize"),drop=F],permu=999)
  sink(file=paste(opt$out,"e_Envfitting2PCA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e05_EnvPCA.pdf",sep=''),width=8, height=5)
  plot(pca,display="sites",xlab=paste("CA1:",round(pca$CA$eig[1]/sum(pca$CA$eig),2)),
       ylab=paste("CA2:",round(pca$CA$eig[2]/sum(pca$CA$eig),2)))
  plot(ef,p.max=0.05)
  dev.off()
  
  pca.c <- rda(dfc~group+Condition(Gender),c,scale=T)
  plot(pca.c,display="sites",xlab=paste("PCA1 ",round(pca.c$CA$eig[1]/sum(pca.c$CA$eig),2)),
       ylab=paste("PCA2 ",round(pca.c$CA$eig[2]/sum(pca.c$CA$eig),2)))
  with(c,ordiellipse(pca.c,group,kind="se",conf=0.95))
  with(c,ordispider(pca.c,group,col="blue",label="T"))
  with(c,ordihull(pca.c,group,col="blue",lty=2))
  
  ca <- cca(dfc)
  pdf(paste(opt$out,"e06_CA.pdf",sep=''),width=8, height=5)
  plot(ca)
  #text(ca,display="sites")
  dev.off()
  ef <- envfit(ca,q[,-which(colnames(q)=="libsize"),drop=F],permu=999)
  sink(file=paste(opt$out,"e_Envfitting2CA.txt",sep='')) 
  ef 
  sink(NULL) 
  cols <- c("blue", "orange", "green","red")
  pdf(paste(opt$out,"e07_EnvCA.pdf",sep=''),width=8, height=5)
  v1 <- round(ca$CA$eig[1]/sum(ca$CA$eig),2)
  v2 <- round(ca$CA$eig[2]/sum(ca$CA$eig),2)
  plot(ca,xlab=paste("CA1:",v1),ylab=paste("CA2:",v2),type="n")
  with(c, points(ca, display = "sites", col = cols[group],pch=c(4,16)[as.numeric(LLI)],bg=cols[group]))
  with(c,ordiellipse(ca,group,kind="se",conf=0.95))
  with(c,ordispider(ca,group,label=TRUE))
  #with(c,ordihull(ca,group,lty=2))
  plot(ef,p.max=0.05)
  dev.off()
  
  
  ca.c <- cca(dfc~group+Condition(Gender),c)
  ef <- envfit(ca.c,q[,-which(colnames(q)=="libsize"),drop=F],permu=999,na.rm=T)
  cols <- c("blue", "orange", "green","red")
  pdf(paste(opt$out,"e07_EnvCA_C.pdf",sep=''),width=8, height=5)
  v1 <- ''#round(summary(ca.c)$concont$importance[2,1],2)
  v2 <- ''#round(summary(ca.c)$concont$importance[2,2],2)
  plot(ca.c,xlab=paste("CA1 ",v1),ylab=paste("CA2 ",v2),type="n")
  with(c, points(ca.c, display = "sites", col = cols[group],pch=c(4,16)[as.numeric(LLI)],bg=cols[group]))
  with(c,ordiellipse(ca.c,group,kind="se",conf=0.95))
  with(c,ordispider(ca.c,group,label=TRUE))
  #with(c,ordihull(ca.c,group,lty=2))
  plot(ef,p.max=0.05)
  dev.off()
  
  dca <- decorana(dfc)
  pdf(paste(opt$out,"e08_DCA.pdf",sep=''),width=8, height=5)
  plot(dca)
  #text(dca,display="sites")
  dev.off()
  ef <- envfit(dca,q[,-which(colnames(q)=="libsize"),drop=F],permu=999)
  sink(file=paste(opt$out,"e_Envfitting2DCA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e09_EnvDCA.pdf",sep=''),width=8, height=5)
  plot(dca,display="sites")
  plot(ef,p.max=0.05)
  #with(c,ordiellipse(dca,group,kind="se",conf=0.95))
  #with(c,ordispider(dca,group,col="blue",label="T"))
  #with(c,ordihull(dca,group,col="blue",lty=2))
  dev.off()
  ## factors
  #ef <- envfit(ca,y,permu=999)
  #ef
  #plot(ca,display="sites",type="p")
  #with(y,ordiellipse(ca,Textura,kind="se"))
  #with(y,ordispider(ca,Textura,col="blue",label=T))
  #with(y,ordihull(ca,Textura,col="blue",lty=2))
  #plot(ef)
  ##
  message("Ordination plots successfuly generated!\n")
}
if (opt$model==T){
  mod0 <- cca(dfc ~ 1,q)
  mod1 <- cca(dfc ~ .,q)
  message("Building a model...")
  #modrev <- step(mod1,scope=list(lower=formula(mod0),upper=formula(mod1)),trace=0) 
  sink(file=paste(opt$out,"m01_fittingmodel.txt",sep='')) 
  mod <- step(mod0,scope=formula(mod1),test="perm")
  
  print("\n Variance Inflation Factors ")
  vif.cca(mod)
  
  print("\n Permutation test - significance of model ")
  anova(mod)
  
  print("\n Permutation test - significance of the terms ")
  anova(mod,by="terms",perm=1000)
  sink(NULL) 
  
  terms <- attr(mod$terms,"term.labels")
  y <- vif.cca(mod)
  Rterms <- c()
  message("Builded-model correction...\n")
  while(length(terms)>0){
    Rterms <- c(Rterms,names(which(y==max(y))))
    terms <- terms[-which(y==max(y))]
    if(length(terms)!=0){
      mody <- cca(as.formula(paste("dfc~",paste(terms,collapse='+'))),q)
      y <- vif.cca(mody)
    }
    modx <- cca(as.formula(paste("dfc~",paste(Rterms,collapse='+'))),q)
    a <- anova(modx,by="terms",perm=1000)
    if (a$`Pr(>F)`[length(a$`Pr(>F)`)-1]>0.1){Rterms <- Rterms[-length(Rterms)]}
  }
  modc <- cca(as.formula(paste("dfc~",paste(Rterms,collapse='+'))),q)
  sink(file=paste(opt$out,"m02_fittedmodel_corrected.txt",sep=''))
  
  print("\n Variance Inflation Factors ")
  vif.cca(modc)
  
  print("\n Permutation test - significance of model ")
  anova(modc)
  
  print("\n Permutation test - significance of the terms ")
  anova(modc,by="terms",perm=1000)

  print("\n Permutation test - marginal effects ")
  anova(modc,by="margin",perm=500)
  
  print("\n Permutation test - significance of the axis ")
  anova(modc,by="axis",perm=1000)
  sink(NULL)
  
  pdf(paste(opt$out,"m01_CCA.pdf",sep=''),width=8, height=5)
  plot(modc)
  text(modc,display="sites")
  dev.off()
  ef <- envfit(ca,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2CA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"m02_lcCCA.pdf",sep=''),width=8, height=5)
  plot(modc,display=c("bp","lc"),xlab=paste("CA1:",round(modc$CCA$eig[1]/sum(modc$CCA$eig),2)),
       ylab=paste("CA2:",round(modc$CCA$eig[2]/sum(modc$CCA$eig),2)))
  dev.off()
  
  spenvcor(modc)
  pdf(paste(opt$out,"m03_waCCA.pdf",sep=''),width=8, height=5)
  plot(modc,display=c("bp","wa"),xlab=paste("CA1:",round(modc$CCA$eig[1]/sum(modc$CCA$eig),2)),
       ylab=paste("CA2:",round(modc$CCA$eig[2]/sum(modc$CCA$eig),2)))
  dev.off()
  
  plot(modc,display="wa")
  #plot(modc,display=c("bp","lc","wa"),type="points")
  #plot(modc,display=c("wa"),type="points")
  #plot(modc,display=c("bp"),type="points")
  #plot(modc,display=c("lc","wa"),type="points")
  #plot(modc,display=c("bp","lc"),type="points")
  #plot(modc,display=c("bp","wa"))
  #ordispider(modc,col="blue")
  
  #pred <- calibrate(modc)
  #head(pred)
  #with(q,plot(Salinity,pred[,"Salinity"]-Salinity,ylab="Prediction Error"))  
  #abline(h=0,col="grey")
  
  
}

## UNDER CONSTRUCTION ## 
#  if(opt$assess_model==''){
#    
#}
### constraint analysis
#dfc ~ Salinity + CT + Arcilla + CO + CE + Cu + Limo + pH
#cca <- cca(dfc ~ Salinity+Fe+CO+CT,x)

#cca <- cca(dfc ~ Salinity + CT + Arcilla,x)
#cca
#plot(cca)
#cca <- cca(dfc ~ Textura,y)
#cca
#plot(cca)
 
#rda2 <- rda(dfc ~ Salinity + CT + Arcilla,x)
#rda2
#plot(rda2,display="sites")
#plot(rda2)

## permutation test

#cca <- cca(dfc ~ Salinity+CT+Arena,x)
#cca
#plot(cca)
#anova(cca,by="term",step=200)
#anova(cca,by="margin",perm=500)
#anova(cca,by="axis",perm=1000)

#mod1 <- cca(dfc ~ .,x)
#mod1
#plot(procrustes(cca(dfc),mod1))

#mod2 <- cca(dfc ~ Salinity + CT + Arcilla,x)
#mod2
#plot(procrustes(cca(dfc),mod2))

#mod3 <- cca(dfc ~ Salinity+Condition(Arcilla),x)
#mod3
#anova(mod3,perm.max=2000)
#plot(mod3,display="sites")
#plot(procrustes(cca(dfc),mod3))

#ef <- envfit(ca,x,permu=999)
#ef
#plot(ca,display="sites")
#plot(ef,p.max=0.05)#

#efc <- envfit(mod2,x,permu=999)
#efc
#plot(mod2)
#plot(ca)
#plot(mod2,display="sites")
#plot(efc,p.max=0.05)

#mod0 <- cca(dfc ~ 1,x)
#mod <- step(mod0,scope=formula(mod1),test="perm")

#vif.cca(mod1)
#vif.cca(mod)

#ccam <- cca(dfc~ Salinity+CO+Arcilla+CT+CE+Cu+Limo,x)
#ccam
#plot(ccam)
#anova(ccam,by="term",step=200)
#anova(ccam,by="margin",perm=500)
#anova(ccam,by="axis",perm=1000)
## weights

#spenvcor(mod)
#spenvcor(mod2)

#k <- cca(dfc ~ Textura,y)
#plot(k,display=c("lc","wa"),type="p")
#ordispider(k,col="blue")

#with(y,ordiellipse(ca,Textura,kind="se"))
#with(y,ordispider(ca,Textura,col="blue",label=T))
#with(y,ordihull(ca,Textura,col="blue",lty=2))
#plot(ef)
#

#barplot(cca$CCA$eig,border=NA,col="gray80",names=colnames(cca$CCA$eig))


