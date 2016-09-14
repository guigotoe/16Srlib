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
# Rscript alpha_div.m.R /path/dataF.rds association_variable /path/outfolder/
# Rscript alpha_div.m.R ~/16Srlib_test/results/dataF.rds Salinity_InterstitialWater ~/16Srlib_test/results/
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
packages(c("metagenomeSeq","vegan","ggplot2","RColorBrewer","RAM","PoiClaClu","zCompositions","reshape2"))

###### end ######

#* input *

f <- commandArgs()[6] # paste(p,'results/dataF.rds',sep='') #
vs <- commandArgs()[7]# 'Salinity,Textura' #
vs <- unlist(strsplit(vs,','))
o <- commandArgs()[8] # paste(p,'results/',sep='') #
## ##
df <- readRDS(f)
dfc <- MRcounts(df,norm=T)
dfc.t <- t(MRcounts(df,norm=T))
df.n <- newMRexperiment(dfc,phenoData=AnnotatedDataFrame(pData(df)),featureData=AnnotatedDataFrame(fData(df)))
## proportion calculation
dfc.p <- t(cmultRepl(t(dfc),method="SQ",output="prop"))
dfp <- newMRexperiment(dfc.p,phenoData=AnnotatedDataFrame(pData(df)),featureData=AnnotatedDataFrame(fData(df)))

## diversity calculation: ##
message("Index calculation...")
indexcal <- data.frame(shannon=round(vegan::diversity(dfc.t),2),
                simpson=round(vegan::diversity(dfc.t, "simpson"),2),
                invsimpson=invsimp <- round(vegan::diversity(dfc.t, "inv"),2))
indexcal$id <- rownames(indexcal)
index <- merge(pData(df),indexcal,by.x=colnames(pData(df))[1],by.y="id")
Ren <- renyi(dfc.t)
rankabund <- radfit(round(dfc.t))

write.table(index,file=paste(o,"index.txt",sep=''),sep="\t",quote=F,row.names=F)
write.table(round(Ren,2),file=paste(o,"renyi.txt",sep=''),sep="\t",quote=F,row.names=F)
sink(file=paste(o,"rankabund.txt",sep='')) 
summary(rankabund) 
sink(NULL) 

## ** PLOTS ** ##
message("Ploting...")

pdf(paste(o,"renyi_plot.pdf",sep=''),width=8, height=5)
plot(Ren,main="Renyi diversities")
dev.off()

pdf(paste(o,"rankabund_plot.pdf",sep=''),width=8, height=5)
plot(rankabund,main="Ranked abundance distribution models",pch=20)
dev.off()

colourCount =NROW(index)
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 9) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount))[-c(1)])
#pie(rep(1,38), col=palette)
pdf(paste(o,"rarefaction_curve.pdf",sep=''),width=8, height=5)
rare <- rarecurve(round(dfc.t), step = 100,col=palette,lwd=3,label=F,ylab="OTUs",xlab="Reads sampled",main="Rarefaction curve")
legend(x=5, y=4,rownames(dfc.t), pch=21,col="#777777", pt.bg=palette,pt.cex=2,cex=.8,bty="n",ncol=9,xjust=0,x.intersp=1,text.width=1515)
dev.off()

for (v in vs){
  shplot <- ggplot(index,aes(x=factor(index[[v]]),y=shannon)) + 
    geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
    ggtitle(paste("Alpha diversity - Shannon index",sep=""))+ylab("Alpha diversity (Shannon entropy)")+
    xlab(v)+scale_fill_discrete(guide=FALSE)+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          plot.title = element_text(lineheight=.12, face="bold"))
  ggsave(paste('shannon_',v,'.pdf'),plot=shplot,path=o,width=8,height=5,device="pdf")
  
  simplot <- ggplot(index,aes(x=factor(index[[v]]),y=simpson)) + 
    geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
    ggtitle(paste("Alpha diversity - Simpson index",sep=""))+ylab("Alpha diversity (Simpson)")+
    xlab(v)+scale_fill_discrete(guide=FALSE)+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          plot.title = element_text(lineheight=.12, face="bold"))
  ggsave(paste('simpson_',v,'.pdf'),plot=simplot,path=o,width=8,height=5,device="pdf")
  
  invsplot <- ggplot(index,aes(x=factor(index[[v]]),y=invsimpson)) + 
    geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
    ggtitle(paste("Alpha diversity - inverse Simpson index",sep=""))+ylab("Alpha diversity (inverse Simpson)")+
    xlab(v)+scale_fill_discrete(guide=FALSE)+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          plot.title = element_text(lineheight=.12, face="bold"))
  ggsave(paste('invsimpson_',v,'.pdf'),plot=invsplot,path=o,width=8,height=5,device="pdf")


## OTU abundance analysis ##
# Compositional Data - proportions

  taxprop(dfp,v,'Phylum',o)
  taxprop(dfp,v,'Class',o)
  taxprop(dfp,v,'Order',o)
  taxprop(dfp,v,'Order',o,u=T)
  taxprop(dfp,v,'Family',o)
  taxprop(dfp,v,'Family',o,u=T)
  taxprop(dfp,v,'Genus',o)
  taxprop(dfp,v,'Genus',o,u=T)
  cl=factor(pData(df.n)[[v]])
  colrs <- colorRampPalette(brewer.pal(8, "Set1"))(9)[as.integer(factor(cl))]
  pdf(paste(o,"MDS_",v,".pdf",sep=''),width=8, height=5)
  plotOrd(df.n,tran=T,usePCA=F,useDist=T,bg=colrs,pch=21)
  legend("bottomleft",levels(cl),box.col=NA,text.col=colorRampPalette(brewer.pal(8, "Set1"))(9),
         title=v,title.col='black')
  dev.off()
}

### Structural analysis ##
# replicates correlation ##
# calculate the pearson correlation between samples in order to see replication issues.

eucldist <- dist(t(dfc))
dismatplot(eucldist,pData(df)[[v]],"Euclidian distance matrix",o)
poisd <- PoissonDistance(t(dfc))
dismatplot(poisd$dd,pData(df)[[v]],"Pearson correlation matrix",o)
brayCurtis <- vegdist(t(dfc))
dismatplot(brayCurtis,pData(df)[[v]],"Bray-Curtis dissimilarity matrix",o)

###
source(toolbox)
#v <- vs[1]
#trials <- pData(df.n)[[v]]

#heatmapColColors <- RAM.pal(cols.needed=length(levels(as.factor(trials))))[as.integer(factor(trials))]
#heatmapCols <- colorRampPalette(brewer.pal(9,"RdBu"))(50)
#plotMRheatmap(obj=df.n,n=200,cexRow=0.4,cexCol=1,trace="none",
#              col=heatmapCols,ColSideColors=heatmapColColors)
#plotCorr(obj=df.n,n=100,cexRow=0.25,cexCol=0.25,trace="none",
#              col=heatmapCols,dendogram="none")
##
#res <- plotRare(df.n,cl=cl,pch=21,bg=cl)
#tmp <- lapply(levels(cl),function(lv) lm(res[,"ident"]~res[,"libSize"] - 1,subset=cl==lv))
#for(i in 1:length(levels(cl))){
#  abline(tmp[[i]],col=i)
#}
#legend("topleft",levels(cl),text.col=c(1,2,3,4),box.col=NA)
##

message("\n ** Index files and plots were succesufy created **\n")
