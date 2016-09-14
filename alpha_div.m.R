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
toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)
#p <- '/home/torres/ikmb_storage/projects/16Srlib_test/'
p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'
packages(c("metagenomeSeq","vegan","ggplot2","RColorBrewer","RAM","PoiClaClu","zCompositions","reshape2"))

###### end ######

#* input *

f <- paste(p,'results/dataF.rds',sep='')#commandArgs()[6] #'/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/results/dataF.rds' # commandArgs()[6] #
v <- 'Salinity' #commandArgs()[7]
o <- paste(p,'results/',sep='') #commandArgs()[8]
## ##
df <- readRDS(f)
dfc <- MRcounts(df,norm=T)
dfc.t <- t(MRcounts(df,norm=T))

df.n <- newMRexperiment(dfc,phenoData=AnnotatedDataFrame(pData(df)),featureData=AnnotatedDataFrame(fData(df)))

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

shplot <- ggplot(index,aes(x=factor(index[[v]]),y=shannon)) + 
  geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - Shannon index",sep=""))+ylab("Alpha diversity (Shannon entropy)")+
  xlab(v)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave(paste('shannon.',v,'.pdf'),plot=shplot,path=o,width=10,height=8,device="pdf")

simplot <- ggplot(index,aes(x=factor(index[[v]]),y=simpson)) + 
  geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - Simpson index",sep=""))+ylab("Alpha diversity (Simpson)")+
  xlab(v)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave(paste('simpson.',v,'.pdf'),plot=simplot,path=o,width=10,height=8,device="pdf")

invsplot <- ggplot(index,aes(x=factor(index[[v]]),y=invsimpson)) + 
  geom_boxplot(aes(fill=factor(index[[v]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - inverse Simpson index",sep=""))+ylab("Alpha diversity (inverse Simpson)")+
  xlab(v)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave(paste('invsimpson.',v,'.pdf'),plot=invsplot,path=o,width=10,height=8,device="pdf")

pdf(paste(o,"renyi_plot.pdf",sep=''),width=10, height=8)
plot(Ren,main="Renyi diversities")
dev.off()

pdf(paste(o,"rankabund_plot.pdf",sep=''),width=10, height=8)
plot(rankabund,main="Ranked abundance distribution models",pch=20)
dev.off()

colourCount =NROW(index)
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 9) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount))[-c(1)])
#pie(rep(1,38), col=palette)
pdf(paste(o,"rarefaction_curve.pdf",sep=''),width=10, height=8)
rare <- rarecurve(round(dfc.t), step = 100,col=palette,lwd=3,label=F,ylab="OTUs",xlab="Reads sampled",main="Rarefaction curve")
legend(x=5, y=4,rownames(dfc.t), pch=21,col="#777777", pt.bg=palette,pt.cex=2,cex=.8,bty="n",ncol=9,xjust=0,x.intersp=1,text.width=1515)
dev.off()
## end ##

### Structural analysis ##
#v <- 'Gender'

## replicates correlation ##
# calculate the pearson correlation between samples in order to see replication issues.

eucldist <- dist(t(dfc))
dismatplot(eucldist,pData(df.n)[[v]],"Euclidian distance matrix")
poisd <- PoissonDistance(t(dfc))
dismatplot(poisd$dd,pData(df.n)[[v]],"Pearson correlation matrix")
brayCurtis <- vegdist(t(dfc))
dismatplot(brayCurtis,pData(df.n)[[v]],"Bray-Curtis dissimilarity matrix")

## OTU abundance analysis ##

vs <- list("Salinity","Textura")
g <- as.data.frame(dfc)
g$id <- rownames(dfc)
gx <- merge(g,fData(df.n),by.x="id",by.y="OTU")
# Compositional Data - proportions
dfc.p <- t(cmultRepl(t(dfc),method="SQ",output="prop"))



dfx <- subset(dfc.p, !is.na(gx$Phylum))
gxx <- subset(gx, !is.na(gx$Phylum))
rownames(dfx) <- gxx$Class
dfc.pt <- as.data.frame(t(dfx))
rownames(dfc.p) <- gx$Phylum
dfc.pt <- as.data.frame(t(dfc.p))
dfc.pt[,is.na(colnames(dfc.pt))==TRUE] <- NULL
dfc.pta <- as.data.frame(sapply(unique(names(dfc.pt)[duplicated(names(dfc.pt))]), function(x) Reduce("+", dfc.pt[ , grep(x, names(dfc.pt))])))
gm <- cbind(dfc.pta,pData(df.n)[vs])
gm_m <- melt(gm,id.vars=vs)
gmx1 <- aggregate(as.formula(paste("value~",'variable+',paste(vs,collapse="+"))),data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
limit=0.01
taxa <- unlist(apply(gmx1,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
gmx1$"taxa" <- as.factor(taxa)
gmx2 <- aggregate(as.formula(paste("value~",'variable+taxa+',paste(vs,collapse="+"))),data=gmx1,FUN=sum) # Because some taxa now are xOthers so we need to summ their values

## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
taxorder <- data.frame(taxa=levels(gmx1$taxa),rate=rep(0,length(levels(gmx1$taxa))))
for (i in levels(gmx2$taxa)){
  a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
  b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
  taxorder$rate[taxorder$taxa==i] <- sum(a,b)
}
taxorder <- taxorder[order(-taxorder[,"rate"]),]
taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
gmx2 <- gmx2[order(gmx2[[v]],gmx2$taxa),]

colourCount = length(levels(gmx2$taxa))
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 3) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount-2))[-c(1)])

taxlevel <- "Phyllum"
title <- "Phyllum abundance distribution"

ggplot(gmx2,aes(x=as.factor(Salinity),y=value,fill=taxa))+geom_boxplot()+
  scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
  scale_x_discrete("Salinity")+ylab("Proportion")+
  ggtitle(paste(title,sep=""))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size=10),
        legend.title = element_text(size=16, face="bold"),
        plot.title = element_text(lineheight=.12, face="bold"))
###
###
dfx <- subset(dfc.p, gx$Phylum!='unclassified' & !is.na(gx$Phylum))
gxx <- subset(gx, gx$Phylum!='unclassified' & !is.na(gx$Phylum))
rownames(dfx) <- gxx$Class
dfc.pt <- as.data.frame(t(dfx))
dfc.pt[,is.na(colnames(dfc.pt))==TRUE] <- NULL
dfc.pta <- as.data.frame(sapply(unique(names(dfc.pt)[duplicated(names(dfc.pt))]), function(x) Reduce("+", dfc.pt[ , grep(x, names(dfc.pt))])))
gm <- cbind(dfc.pta,pData(df.n)[vs])
gm_m <- melt(gm,id.vars=vs)
gmx1 <- aggregate(value ~ Punto+variable+Salinity,data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
limit=0.01
taxa <- unlist(apply(gmx1,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
gmx1$"taxa" <- as.factor(taxa)
gmx2 <- aggregate(value ~ Punto+taxa+Salinity,data=gmx1,FUN=sum) # Because some taxa now are xOthers so we need to summ their values

## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
taxorder <- data.frame(taxa=levels(gmx1$taxa),rate=rep(0,length(levels(gmx1$taxa))))
for (i in levels(gmx2$taxa)){
  a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
  b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
  taxorder$rate[taxorder$taxa==i] <- sum(a,b)
}
taxorder <- taxorder[order(-taxorder[,"rate"]),]
taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
gmx2 <- gmx2[order(gmx2[[v]],gmx2$taxa),]

colourCount = length(levels(gmx2$taxa))
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 3) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount-2))[-c(1)])

taxlevel <- "Class"
title <- "Class abundance distribution"

ggplot(gmx2,aes(x=as.factor(Salinity),y=value,fill=taxa))+geom_boxplot()+
  scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
  scale_x_discrete("Salinity")+ylab("Proportion")+
  ggtitle(paste(title,sep=""))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size=10),
        legend.title = element_text(size=16, face="bold"),
        plot.title = element_text(lineheight=.12, face="bold"))


###
###
dfx <- subset(dfc.p, gx$Phylum!='unclassified'& gx$Order!='unclassified' & !is.na(gx$Phylum))
gxx <- subset(gx, gx$Phylum!='unclassified' & gx$Order!='unclassified' & !is.na(gx$Phylum))
rownames(dfx) <- gxx$Order
dfc.pt <- as.data.frame(t(dfx))
dfc.pt[,is.na(colnames(dfc.pt))==TRUE] <- NULL
dfc.pta <- as.data.frame(sapply(unique(names(dfc.pt)[duplicated(names(dfc.pt))]), function(x) Reduce("+", dfc.pt[ , grep(x, names(dfc.pt))])))
gm <- cbind(dfc.pta,pData(df.n)[vs])
gm_m <- melt(gm,id.vars=vs)
gmx1 <- aggregate(value ~ Punto+variable+Salinity,data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
limit=0.01
taxa <- unlist(apply(gmx1,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
gmx1$"taxa" <- as.factor(taxa)
gmx2 <- aggregate(value ~ Punto+taxa+Salinity,data=gmx1,FUN=sum) # Because some taxa now are xOthers so we need to summ their values

## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
taxorder <- data.frame(taxa=levels(gmx1$taxa),rate=rep(0,length(levels(gmx1$taxa))))
for (i in levels(gmx2$taxa)){
  a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
  b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
  taxorder$rate[taxorder$taxa==i] <- sum(a,b)
}
taxorder <- taxorder[order(-taxorder[,"rate"]),]
taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
gmx2 <- gmx2[order(gmx2[[v]],gmx2$taxa),]

colourCount = length(levels(gmx2$taxa))
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 3) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount-2))[-c(1)])

taxlevel <- "Order"
title <- "Order abundance distribution"

ggplot(gmx2,aes(x=as.factor(Salinity),y=value,fill=taxa))+geom_boxplot()+
  scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
  scale_x_discrete("Salinity")+ylab("Proportion")+
  ggtitle(paste(title,sep=""))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size=10),
        legend.title = element_text(size=16, face="bold"),
        plot.title = element_text(lineheight=.12, face="bold"))
###
##

dfx <- subset(dfc.p, gx$Phylum!='unclassified'& gx$Order!='unclassified' & !is.na(gx$Phylum))
gxx <- subset(gx, gx$Phylum!='unclassified' & gx$Order!='unclassified' & !is.na(gx$Phylum))
rownames(dfx) <- gxx$Genus
dfc.pt <- as.data.frame(t(dfx))
dfc.pt[,is.na(colnames(dfc.pt))==TRUE] <- NULL
dfc.pta <- as.data.frame(sapply(unique(names(dfc.pt)[duplicated(names(dfc.pt))]), function(x) Reduce("+", dfc.pt[ , grep(x, names(dfc.pt))])))
gm <- cbind(dfc.pta,pData(df.n)[vs])
gm_m <- melt(gm,id.vars=vs)
gmx1 <- aggregate(value ~ Punto+variable+Salinity,data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
limit=0.01
taxa <- unlist(apply(gmx1,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
gmx1$"taxa" <- as.factor(taxa)
gmx2 <- aggregate(value ~ Punto+taxa+Salinity,data=gmx1,FUN=sum) # Because some taxa now are xOthers so we need to summ their values

## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
taxorder <- data.frame(taxa=levels(gmx1$taxa),rate=rep(0,length(levels(gmx1$taxa))))
for (i in levels(gmx2$taxa)){
  a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
  b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
  taxorder$rate[taxorder$taxa==i] <- sum(a,b)
}
taxorder <- taxorder[order(-taxorder[,"rate"]),]
taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
gmx2 <- gmx2[order(gmx2[[v]],gmx2$taxa),]

colourCount = length(levels(gmx2$taxa))
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 3) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount-2))[-c(1)])

taxlevel <- "Genus"
title <- "Genus abundance distribution"

ggplot(gmx2,aes(x=as.factor(Salinity),y=value,fill=taxa))+geom_boxplot()+
  scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
  scale_x_discrete("Salinity")+ylab("Proportion")+
  ggtitle(paste(title,sep=""))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size=10),
        legend.title = element_text(size=16, face="bold"),
        plot.title = element_text(lineheight=.12, face="bold"))

###

dfx <- subset(dfc.p, gx$Phylum!='unclassified'& gx$Family!='unclassified' & !is.na(gx$Phylum))
gxx <- subset(gx, gx$Phylum!='unclassified' & gx$Family!='unclassified' & !is.na(gx$Phylum))
rownames(dfx) <- gxx$Family
dfc.pt <- as.data.frame(t(dfx))
dfc.pt[,is.na(colnames(dfc.pt))==TRUE] <- NULL
dfc.pta <- as.data.frame(sapply(unique(names(dfc.pt)[duplicated(names(dfc.pt))]), function(x) Reduce("+", dfc.pt[ , grep(x, names(dfc.pt))])))
gm <- cbind(dfc.pta,pData(df.n)[vs])
gm_m <- melt(gm,id.vars=vs)
gmx1 <- aggregate(value ~ Punto+variable+Salinity,data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
limit=0.01
taxa <- unlist(apply(gmx1,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
gmx1$"taxa" <- as.factor(taxa)
gmx2 <- aggregate(value ~ Punto+taxa+Salinity,data=gmx1,FUN=sum) # Because some taxa now are xOthers so we need to summ their values

## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
taxorder <- data.frame(taxa=levels(gmx1$taxa),rate=rep(0,length(levels(gmx1$taxa))))
for (i in levels(gmx2$taxa)){
  a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
  b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
  taxorder$rate[taxorder$taxa==i] <- sum(a,b)
}
taxorder <- taxorder[order(-taxorder[,"rate"]),]
taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
gmx2 <- gmx2[order(gmx2[[v]],gmx2$taxa),]

colourCount = length(levels(gmx2$taxa))
base <- colorRampPalette(brewer.pal(8, "Set1"))(9)[1:3]
if (colourCount <= 3) {palette <- base
}else palette <- c(base,RAM.pal(cols.needed=(colourCount-2))[-c(1)])

taxlevel <- "Family"
title <- "Family abundance distribution"

ggplot(gmx2,aes(x=as.factor(Salinity),y=value,fill=taxa))+geom_boxplot()+
  scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
  scale_x_discrete("Salinity")+ylab("Proportion")+
  ggtitle(paste(title,sep=""))+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position="bottom",legend.box="horizontal",
        legend.text = element_text(size=10),
        legend.title = element_text(size=16, face="bold"),
        plot.title = element_text(lineheight=.12, face="bold"))

#if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,idslab,'.pdf',sep=""),width=12, height=8)

trials <- pData(df.n)[[v]]

heatmapColColors <- RAM.pal(cols.needed=length(levels(as.factor(trials))))[as.integer(factor(trials))]
heatmapCols <- colorRampPalette(brewer.pal(9,"RdBu"))(50)
plotMRheatmap(obj=df.n,n=200,cexRow=0.4,cexCol=1,trace="none",
              col=heatmapCols,ColSideColors=heatmapColColors)
plotCorr(obj=df.n,n=100,cexRow=0.25,cexCol=0.25,trace="none",
              col=heatmapCols,dendogram="none")
#v <- 'Gender'
#x <- cut(pData(df.n)[["Age"]],breaks=9)
cl=factor(pData(df.n)[[v]])
colrs <- colorRampPalette(brewer.pal(8, "Set1"))(9)[as.integer(factor(cl))]
plotOrd(df.n,tran=T,usePCA=F,useDist=T,bg=colrs,pch=21)
legend("bottomright",levels(cl),box.col=NA,text.col=colorRampPalette(brewer.pal(8, "Set1"))(9))



##
res <- plotRare(df.n,cl=cl,pch=21,bg=cl)

tmp <- lapply(levels(cl),function(lv) lm(res[,"ident"]~res[,"libSize"] - 1,subset=cl==lv))
for(i in 1:length(levels(cl))){
  abline(tmp[[i]],col=i)
}
legend("topleft",levels(cl),text.col=c(1,2,3,4),box.col=NA)
##


message("\n ** Index files and plots were succesufy created **\n")



tax_graph <- function(df,limit=0.06,taxlevel="Genus",savef=NULL,ids=F,method="SQ",taxon=NULL){
  gm <- fillzeros(df,output="prop",method=method)
  if("Bacteroidetes"%in%colnames(gm)){gm$"taxlab" <- gm[,"Firmicutes"]  # phyllum reference to graph order
  }else gm$"taxlab" <- rep(1,NROW(gm))
  colnames(gm)[which(names(gm) == "gut")] <- "RC9_gut_group"
  colnames(gm)[which(names(gm) == "Incertae")] <- "Incertae_Sedis"
  gm_m <- melt(gm,id.vars=c("id","gender","age","agex","age.group","age.mean","taxlab"))
  ranks <- unlist(lapply(gm_m$"age", function(x) getrank(x)))
  taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
  gm_m$"ranks" <- as.factor(ranks)
  gm_m$"taxa" <- as.factor(taxa) 
  gmx <- within(gm_m,position <- factor(age,levels=names(sort(table(age)))))
  gmx2 <- aggregate(value ~ id+taxa+age+agex+position+ranks+age.group+age.mean+gender+taxlab,data=gmx,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
  
  ## before plot we need to order the taxa levels according their values of abundance and prevalence in individuals
  taxorder <- data.frame(taxa=levels(gmx2$taxa),rate=rep(0,length(levels(gmx2$taxa))))
  for (i in levels(gmx2$taxa)){
    a <- sum(unlist(lapply(gmx2$value[gmx2$taxa==i],sum))) # OTU global abundace. max = No. of individuals 
    b <- table(gmx2$taxa)[i][[1]]                          # OTU prevalence. max = No. of individuals
    taxorder$rate[taxorder$taxa==i] <- sum(a,b)
  }
  taxorder <- taxorder[order(-taxorder[,"rate"]),]
  taxorder <- taxorder[-which(taxorder$taxa=="Others"),]
  gmx2$taxa <- factor(gmx2$taxa,levels=c(as.character(taxorder$taxa),"Others"))
  gmx2 <- gmx2[order(gmx2$agex,gmx2$taxa),]
  
  ## placing the colors to the taxas
  colourCount = length(levels(gmx2$taxa))
  base <- colorRampPalette(brewer.pal(8, "Set1"))(colourCount-1)
  base <- base[c(1,3,2,c(4:length(base)))]
  base2 <- colorRampPalette(brewer.pal(8, "Set1"))(9)
  base2 <- base[c(1,3,2,c(4:length(base)))] 
  if (colourCount <= 10) {palette <- c(base,"#3D3D3D")
  }else if (colourCount <= 19) {palette <- c(base2,colorRampPalette(brewer.pal(8, "Accent"))(colourCount-10),"#3D3D3D")
  }else palette <- c(base2,colorRampPalette(brewer.pal(8, "Set3"))(9),
                     colorRampPalette(brewer.pal(8, "Accent"))(colourCount-19),"#3D3D3D")
  
  if(taxlevel=="Phylla"){
    gm$ratioFB <- unlist(apply(gm,1,function(x) round(as.numeric(x["Firmicutes"])/as.numeric(x["Bacteroidetes"]),2)))
    print(aggregate(gm[,"ratioFB"], list(gm$age.group), mean))
    print(c(summary(as.numeric(gm[,"Firmicutes"])),summary(as.numeric(gm[,"Bacteroidetes"])),summary(as.numeric(gm[,"Proteobacteria"]))))
  }
  
  plot <- function(gm.g,title,taxlevel,palette,savef=NULL,ids=F){
    x_labels <- gm[,c("agex","age","id","taxlab")][order(gm[,c("agex","age","id","taxlab")][,"agex"]),]
    if (ids==T){ 
      x_lab <- unlist(apply(x_labels,1,function(x) paste(x[3],x[2],sep="_")))
      idslab <- '_ids'
    }else {
      x_lab <- x_labels$age
      idslab <- ''
    }
    #gmx2.g <- subset(gmx2,gmx2$gender==g) ## each gender
    
    # to be color consistente 
    e <- length(table(gm.g$taxa)[table(gm.g$taxa)==0])
    if (e!=0) {palette.g <- palette[-c((length(palette)-e):(length(palette)-1))]
    }else palette.g <- palette
    
    ggplot(gm.g, aes(x=as.factor(agex),y=value,fill=taxa))+
      geom_bar(stat="identity",width=1)+
      scale_fill_manual(name=taxlevel,values=palette.g) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
      scale_x_discrete("Age",labels=as.character(x_lab))+ylab("Proportion")+
      guides(fill=guide_legend(ncol=6,keywidth=1, keyheight=1))+
      ggtitle(title)+#facet_wrap(~ age+facet_wrap(~ taxa))+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
            panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
            legend.position="bottom",legend.box="horizontal",
            legend.text = element_text(size=14),
            legend.title = element_text(size=16, face="bold"),
            plot.title = element_text(lineheight=.12, face="bold"))
    if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,idslab,'.pdf',sep=""),width=12, height=8)
    
    ggplot(gm.g,aes(x=as.factor(age.group),y=value,fill=taxa))+geom_boxplot()+
      scale_fill_manual(name=taxlevel,values=palette.g) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
      scale_x_discrete("")+ylab("Proportion")+
      ggtitle(paste(title,sep=""))+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
            panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
            legend.position="bottom",legend.box="horizontal",
            legend.text = element_text(size=14),
            legend.title = element_text(size=16, face="bold"),
            plot.title = element_text(lineheight=.12, face="bold"))
    if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,idslab,'_boxplot.pdf',sep=""),width=12, height=8)
    
    if(title=="All" & taxlevel=="Phylla"){
      phylabs <- x_labels[order(x_labels[,"taxlab"]),]
      x_lab <- phylabs$age
      
      # for comparisson effects with paper Claesson 2011
      #x <- levels(as.factor(gm.g$taxlab))
      #y <- levels(as.factor(gm.g$taxa))
      #gm.g$taxlab <- factor(gm.g$taxlab,levels=sort(x,decreasing = T))    
      #gm.g$taxa <- factor(gm.g$taxa,levels=y[c(2,1,c(3:length(y)))])
      #gm.g$taxa <- factor(gm.g$taxa,levels=y)
      gm.g$taxlab <- factor(gm.g$taxlab,levels=sort(x))
      ###
      
      ggplot(gm.g, aes(x=as.factor(taxlab),y=value,fill=taxa))+
        geom_bar(stat="identity",width=1)+
        scale_fill_manual(name=taxlevel,values=palette) + #limits=c(levels(gmx2[1]),levels(gmx2[length(levels(gmx2))]))
        scale_x_discrete("Age",labels=as.character(x_lab))+ylab("Proportion")+
        guides(fill=guide_legend(ncol=6,keywidth=1, keyheight=1))+
        ggtitle(title)+#facet_wrap(~ age+facet_wrap(~ taxa))+
        theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
              panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
              legend.position="bottom",legend.box="horizontal",
              legend.text = element_text(size=14),
              legend.title = element_text(size=16, face="bold"),
              plot.title = element_text(lineheight=.12, face="bold"))
      if(!is.null(savef)) ggsave(paste(savef,title,'_',taxlevel,idslab,'_boxplot_2ClaessonComp.pdf',sep=""),width=12, height=8)
    }
  }
  ## Ploting by genders
  # labels for x axis
  
  for (g in levels(as.factor(gm$gender))){
    title <- paste(toupper(substring(as.character(g),1,1)),substring(as.character(g),2),sep="")
    #gm.g <- subset(gm,gm$gender==g) ## each gender
    gmx2.g <- subset(gmx2,gmx2$gender==g)
    plot(gmx2.g,title,taxlevel=taxlevel,palette=palette,savef=savef)
  }
  plot(gmx2,title="All",taxlevel=taxlevel,palette=palette,savef=savef)
}

