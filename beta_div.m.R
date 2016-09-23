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
# Rscript beta_div.R -- help
# Rscript beta_div.R -i ~/16Srlib_test/results/dataF.rds -o ~/16Srlib_test/results/ -e [options] 
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
packages(c("metagenomeSeq","reshape2","vegan","ggplot2","optparse"))

## Options ##
p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'

option_list <- list(
  make_option(c("-i","--data"),type="character",default=paste(p,'results/dataF.rds',sep=''),
              help="Path to input rds file"),
  make_option(c("-o","--out"),type="character",default=paste(p,'results/',sep=''),
             help="Path to output directory [default %default]"),
  make_option(c("-e","--exploratory"),action="store_true",default=FALSE,
              help="Perform exploratory analysis"),
  make_option(c("-m","--model"),action="store_true",default=FALSE,
              help="Build constrained model based on AIC selection criterion"),
  make_option(c("-am","--assess_model"),type="character",default='',#NULL,
              help="Model's terms assessed by permutation tests; for new model: vars,separated,by,comma"),
  make_option(c("-C","--constraints"),type="character",default=NULL,
              help="Set of constraints used by -coa: vars,separated,by,comma"),
  make_option(c("-coa","--constrained_analysis"),action="store_true",default=FALSE,
              help="Perform constrained ordination analysis using -C constraints"),
  make_option(c("-B","--factors"),type="character",default=NULL,
              help="Set of factors used by -b: vars,separated,by,comma"),
  make_option(c("-b","--beta"),action="store_true",default=FALSE,
              help="Perform beta diversty between -B variable classes"),
  make_option(c("-f","--filter"),type="double",default=0.3,
              help="Percentile as threshold of low abundant features ")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.null(opt$data)){stop(sprintf("There is not file specified"))}

#### Preparing the input data ####
data <- readRDS(opt$data)
dlog <- log2(MRcounts(data)+1)
dlogm <- sort.default(unlist(apply(dlog,1,function(x) mean(x))))
csnorm <- cumsum(dlogm)/sum(dlogm)
totrim <- names(csnorm[which(csnorm<=quantile(csnorm,probs=opt$filter))])
totrim <- match(totrim,rownames(data))
dataTrimed <- data[-totrim,]

dfc0 <- t(MRcounts(dataTrimed,norm=T))
dfc <- t(log2(MRcounts(dataTrimed,norm=T)+1))

q <- pData(dataTrimed)[, !sapply(pData(dataTrimed), is.factor)]
c <- pData(dataTrimed)[, sapply(pData(dataTrimed), is.factor)]
#taxa <- fData(df)[,which(colnames(fData(df))%in%c("Phylum","Class","Order","Family"))] 
#rankindex(scale(x),t(MRcounts(df,norm=T)),c("euc","man","bray","jac","kul"))

###### end ######

if (opt$exploratory==T){
  pdf(paste(opt$out,"e_data_filtering.pdf",sep=''),width=8, height=5)
  plot(csnorm)
  abline(h=quantile(csnorm,probs=opt$filter))
  dev.off()
  #dist <- vegdist(dfc)
  # 1) Direction of the gradient: The arrow points of the direction of most rapid change 
  #    in the environmental variable.
  # 2) Strength of the gradient: The length of the arrow is proportional to the correlation 
  #    between ordination and environmental variable.
  
  nmds <- metaMDS(dfc,trace=F)
  pdf(paste(opt$out,"e_NMDS.pdf",sep=''),width=8, height=5)
  ordiplot(nmds,type="t",display="sites")
  dev.off()
  ef <- envfit(nmds,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2NMDS.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e_EnvNMDS.pdf",sep=''),width=8, height=5)
  plot(nmds,display="sites")
  plot(ef,p.max=0.05)
  dev.off()
  
  pca <- rda(dfc,scale=TRUE)
  pdf(paste(opt$out,"e_PCA.pdf",sep=''),width=8, height=5)
  plot(pca)
  text(pca,display="sites")
  dev.off()
  ef <- envfit(pca,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2PCA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e_EnvPCA.pdf",sep=''),width=8, height=5)
  plot(pca,display="sites",xlab=paste("CA1:",round(pca$CA$eig[1]/sum(pca$CA$eig),2)),
       ylab=paste("CA2:",round(pca$CA$eig[2]/sum(pca$CA$eig),2)))
  plot(ef,p.max=0.05)
  dev.off()
  
  ca <- cca(dfc)
  pdf(paste(opt$out,"e_CA.pdf",sep=''),width=8, height=5)
  plot(ca)
  text(ca,display="sites")
  dev.off()
  ef <- envfit(ca,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2CA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e_EnvCA.pdf",sep=''),width=8, height=5)
  plot(ca,display="sites",xlab=paste("CA1:",round(ca$CA$eig[1]/sum(ca$CA$eig),2)),
       ylab=paste("CA2:",round(ca$CA$eig[2]/sum(ca$CA$eig),2)))
  plot(ef,p.max=0.05)
  dev.off()

  dca <- decorana(dfc)
  pdf(paste(opt$out,"e_CA.pdf",sep=''),width=8, height=5)
  plot(dca)
  text(dca,display="sites")
  dev.off()
  ef <- envfit(dca,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2DCA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"e_EnvDCA.pdf",sep=''),width=8, height=5)
  plot(dca,display="sites")
  plot(ef,p.max=0.05)
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
  sink(file=paste(opt$out,"m_fittingmodel.txt",sep='')) 
  mod <- step(mod0,scope=formula(mod1),test="perm")
  
  print(" Variance Inflation Factors ")
  vif.cca(mod)
  
  print(" Permutation test - significance of model ")
  anova(mod)
  
  print(" Permutation test - significance of the terms ")
  anova(mod,by="terms",perm=1000)
  sink(NULL) 
  
  terms <- attr(mod$terms,"term.labels")
  y <- vif.cca(mod)
  Rterms <- c()
  message("Builded-model correction...")
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
  sink(file=paste(opt$out,"m_fittedmodel_corrected.txt",sep=''))
  
  print(" Variance Inflation Factors ")
  vif.cca(modc)
  
  print(" Permutation test - significance of model ")
  anova(modc)
  
  print(" Permutation test - significance of the terms ")
  anova(modc,by="terms",perm=1000)

  print(" Permutation test - marginal effects ")
  anova(modc,by="margin",perm=500)
  
  print(" Permutation test - significance of the axis ")
  anova(modc,by="axis",perm=1000)
  sink(NULL)
  
  pdf(paste(opt$out,"m_CCA.pdf",sep=''),width=8, height=5)
  plot(modc)
  text(modc,display="sites")
  dev.off()
  ef <- envfit(ca,q,permu=999)
  sink(file=paste(opt$out,"e_Envfitting2CA.txt",sep='')) 
  ef 
  sink(NULL) 
  pdf(paste(opt$out,"m_EnvCCA.pdf",sep=''),width=8, height=5)
  plot(ca,display="sites",xlab=paste("CA1:",round(modc$CCA$eig[1]/sum(modc$CCA$eig),2)),
       ylab=paste("CA2:",round(modc$CCA$eig[2]/sum(modc$CCA$eig),2)))
  plot(modc,display="wa")
  dev.off()
  
  spenvcor(modc)
  plot(modc,display=c("bp","lc","wa"),type="points")
  plot(modc,display=c("lc","wa"),type="points")
  ordispider(modc,col="blue")
  
  pred <- calibrate(modc)
  head(pred)
  with(q,plot(Salinity,pred[,"Salinity"]-Salinity,ylab="Prediction Error"))  
  abline(h=0,col="grey")
  
  plot(mod=)
}
  if(opt$assess_model==''){
    
}
### constraint analysis
dfc ~ Salinity + CT + Arcilla + CO + CE + Cu + Limo + pH
#cca <- cca(dfc ~ Salinity+Fe+CO+CT,x)

cca <- cca(dfc ~ Salinity + CT + Arcilla,x)
cca
plot(cca)


cca <- cca(dfc ~ Textura,y)
cca
plot(cca)
 
rda2 <- rda(dfc ~ Salinity + CT + Arcilla,x)
rda2
plot(rda2,display="sites")
plot(rda2)

## permutation test

cca <- cca(dfc ~ Salinity+CT+Arena,x)
cca
plot(cca)
anova(cca,by="term",step=200)
anova(cca,by="margin",perm=500)
anova(cca,by="axis",perm=1000)

mod1 <- cca(dfc ~ .,x)
mod1
plot(procrustes(cca(dfc),mod1))

mod2 <- cca(dfc ~ Salinity + CT + Arcilla,x)
mod2
plot(procrustes(cca(dfc),mod2))

mod3 <- cca(dfc ~ Salinity+Condition(Arcilla),x)
mod3
anova(mod3,perm.max=2000)
plot(mod3,display="sites")
plot(procrustes(cca(dfc),mod3))

ef <- envfit(ca,x,permu=999)
ef
plot(ca,display="sites")
plot(ef,p.max=0.05)

efc <- envfit(mod2,x,permu=999)
efc
plot(mod2)
plot(ca)
plot(mod2,display="sites")
plot(efc,p.max=0.05)

mod0 <- cca(dfc ~ 1,x)
mod <- step(mod0,scope=formula(mod1),test="perm")

vif.cca(mod1)
vif.cca(mod)

ccam <- cca(dfc~ Salinity+CO+Arcilla+CT+CE+Cu+Limo,x)
ccam
plot(ccam)
anova(ccam,by="term",step=200)
anova(ccam,by="margin",perm=500)
anova(ccam,by="axis",perm=1000)
## weights

spenvcor(mod)
spenvcor(mod2)

k <- cca(dfc ~ Textura,y)
plot(k,display=c("lc","wa"),type="p")
ordispider(k,col="blue")

with(y,ordiellipse(ca,Textura,kind="se"))
with(y,ordispider(ca,Textura,col="blue",label=T))
with(y,ordihull(ca,Textura,col="blue",lty=2))
plot(ef)
#

barplot(cca$CCA$eig,border=NA,col="gray80",names=colnames(cca$CCA$eig))


