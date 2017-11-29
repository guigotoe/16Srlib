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
toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
p <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/2017_11_gg2/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'
packages(c("metagenomeSeq","vegan","ggplot2","RColorBrewer","RAM","PoiClaClu","zCompositions","reshape2"))
###### end ######

#* input *

f <- paste(p,'dataFcp_l_0.2.rds',sep='')#commandArgs()[6] # paste(p,'results/dataF.rds',sep='') #
vs <- 'Sample' #'group,Gender'#commandArgs()[7]#
vs <- unlist(strsplit(vs,','))
extravar <- 'ID_ref'
o <- paste(p,'alpha/',sep='')#commandArgs()[8] # paste(p,'results/',sep='') #
##
if(dir.exists(o)){message('Out-folder already exist, files will be overwritten')
}else dir.create(o,showWarnings=F)

## ##
df <- readRDS(f)
pData(df)$ID_ref <- factor(pData(df)$ID_ref,levels=c('low','med','high'))
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
#index$ID_ref

write.table(index,file=paste(o,"index.txt",sep=''),sep="\t",quote=F,row.names=F)
write.table(round(Ren,2),file=paste(o,"renyi.txt",sep=''),sep="\t",quote=F,row.names=F)
sink(file=paste(o,"rankabund.txt",sep='')) 
summary(rankabund) 
sink(NULL) 

## ** PLOTS ** ##
message("Ploting...")

#pdf(paste(o,"renyi_plot.pdf",sep=''),width=8, height=5,onefile=FALSE)
#plot(Ren,main="Renyi diversities")
#dev.off()

#pdf(paste(o,"rankabund_plot.pdf",sep=''),width=8, height=5,onefile=FALSE)
#plot(rankabund,main="Ranked abundance distribution models",pch=20)
#dev.off()

#colourCount =NROW(index)
#base <- colorRampPalette(brewer.pal(8, "Set1"))(9)
#if (colourCount <= 9) {palette <- base
#}else palette <- c(base,RAM.pal(cols.needed=(colourCount))[-c(1)])
if (length(levels(index[[extravar]]))>9){ l <- length(levels(index[[extravar]]))
}else l <- 9
sampleColors= colorRampPalette(brewer.pal(8, "Set1"))(l)[c(1:length(levels(index[[extravar]])))]
sample_colors_dframe = data.frame(samples=levels(index[[extravar]]), colors=sampleColors)
sample_colors = unlist(lapply(index[[extravar]], function(x) return(as.character(sample_colors_dframe[x,2]))))
#pie(rep(1,10), col=sample_colors)
raredata <- round(dfc.t)
names_raredata <- data.frame(rowids=rownames(raredata),ref=index[match(as.character(index[[1]]),rownames(raredata))][[extravar]])
names_raredata <- names_raredata[with(names_raredata, order(ref)),]
raredata <- raredata[names_raredata$rowids,]
pdf(paste(o,"rarefaction_curve.pdf",sep=''),width=8, height=6,onefile=FALSE)
rare <- rarecurve(raredata,step = 100,col=sample_colors,lwd=3,label=F,ylab="Observed OTU",xlab="Reads sampled",main="Rarefaction curve")
#as.character(names_raredata$ref)
legend("bottomright",as.character(sample_colors_dframe$samples), pch=21,col="#777777",pt.bg=as.character(sample_colors_dframe$colors),pt.cex=1.2,cex=.8,bty="n",ncol=10,xjust=1,x.intersp=0.5,y.intersp=1)#,text.width=1200)
dev.off()
v <- 'ID_ref'
for (v in vs){
  shplot <- ggplot(index,aes(x=factor(index[[v]]),y=shannon)) + 
    geom_boxplot(aes(fill=factor(index[[v]])),width=0.9) +# geom_jitter(position=position_jitter(width=0.3))+
    ggtitle(paste("Alpha diversity - Shannon index",sep=""))+ylab("Alpha diversity (Shannon entropy)")+
    scale_fill_manual(values=as.character(sample_colors_dframe$colors),guide=FALSE)+
    scale_x_discrete(expand=c(0,0.5))+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=20),
          axis.text.y  = element_text(size=20),
          axis.title.y=element_text(size=20,face="bold"),
          axis.title.x=element_blank(),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          plot.title = element_blank(),#element_text(lineheight=.12, face="bold"),
          plot.background= element_blank()
    )
  ggsave(paste('shannon_',v,'.pdf'),plot=shplot,path=o,width=8,height=6,device="pdf")
  
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
  tphylum <- taxprop(dfp,v,'Phylum',o,u=F,z=F,limit=10,legendncol=4,textsize=10)
  tclass <- taxprop(dfp,v,'Class',o,u=F,z=F,legendncol=3,textsize=10)
  torder <- taxprop(dfp,v,'Order',o,u=F,z=F,legendncol=3,textsize=10)
  tfamily <- taxprop(dfp,v,'Family',o,u=F,z=F,legendncol=3,textsize=8)
  tgenus <- taxprop(dfp,v,'Genus',o,u=F,z=F,legendncol=3,textsize=7)
  cl=factor(pData(df.n)[[v]])
  colrs <- colorRampPalette(brewer.pal(8, "Set1"))(9)[as.integer(factor(cl))]
  pdf(paste(o,"MDS_",v,".pdf",sep=''),width=8, height=5,onefile=FALSE)
  plotOrd(df.n,tran=T,usePCA=F,useDist=T,bg=colrs,pch=21)
  legend("bottomleft",levels(cl),box.col=NA,text.col=colorRampPalette(brewer.pal(8, "Set1"))(9),
         title=v,title.col='black')
  dev.off()
}
min(as.numeric(tphylum[tphylum$OTU=='Planctomycetes',1:(ncol(tphylum)-2)]))*100
max(as.numeric(tphylum[tphylum$OTU=='Planctomycetes',1:(ncol(tphylum)-2)]))*100

### Phylum T-test actinobacteria med vs low and med vs high and low vs high (no sig)
otu <- 'Actinobacteria'
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
### Phylum T-test Bacteroidetes med vs low (sig) and med vs high (X) and low vs high (sig)
otu <- 'Bacteroidetes'
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
### Phylum T-test Planctomycetes med vs low (sig) and med vs high (sig) and low vs high (sig)
otu <- 'Planctomycetes'
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
### Phylum T-test Gemmatimonadetes med vs low (x) and med vs high (sig) and low vs high (sig)
otu <- 'Gemmatimonadetes'
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('2A','2B','2C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))
t.test(as.numeric(tphylum[tphylum$OTU==otu,c('3A','3B','3C')]),as.numeric(tphylum[tphylum$OTU==otu,c('4A','4B','4C')]))

### class ###

mean(as.numeric(tclass[tclass$OTU=='Alphaproteobacteria',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Alphaproteobacteria',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Alphaproteobacteria',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Gammaproteobacteria',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Gammaproteobacteria',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Gammaproteobacteria',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Deltaproteobacteria',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Deltaproteobacteria',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Deltaproteobacteria',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Gemmatimonadetes',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Gemmatimonadetes',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Gemmatimonadetes',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Phycisphaerae',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Phycisphaerae',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Phycisphaerae',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Nitrospira',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Nitrospira',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Nitrospira',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='S085',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='S085',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='S085',c('4A','4B','4C')]))*100

mean(as.numeric(tclass[tclass$OTU=='Acidimicrobiia',c('3A','3B','3C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Acidimicrobiia',c('2A','2B','2C')]))*100
mean(as.numeric(tclass[tclass$OTU=='Acidimicrobiia',c('4A','4B','4C')]))*100

### order##

mean(as.numeric(torder[torder$OTU=='Acidimicrobiales',c('3A','3B','3C')]))*100
mean(as.numeric(torder[torder$OTU=='Acidimicrobiales',c('2A','2B','2C')]))*100
mean(as.numeric(torder[torder$OTU=='Acidimicrobiales',c('4A','4B','4C')]))*100

### Structural analysis ##
# replicates correlation ##
# calculate the pearson correlation between samples in order to see replication issues.

#colnames(dfc) <- index[match(colnames(dfc),index$ID),3]
eucldist <- dist(t(dfc))
dismatplot(eucldist,pData(df)[[extravar]],"Euclidian distance matrix",o)
poisd <- PoissonDistance(t(dfc))
dismatplot(poisd$dd,pData(df)[[extravar]],"Pearson correlation matrix",o)
brayCurtis <- vegdist(t(dfc))
dismatplot(brayCurtis,pData(df)[[extravar]],"Bray-Curtis dissimilarity matrix",o)

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
