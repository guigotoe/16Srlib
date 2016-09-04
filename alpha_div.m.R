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

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "toolbox.R")
#source("/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R")
#source("/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R")
source(other.name)

packages(c("metagenomeSeq","vegan","ggplot2","RColorBrewer"))

###### end ######

#* input *

f <- commandArgs()[6] #'/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/results/dataF.rds' # commandArgs()[6] #
assocvar <- commandArgs()[7] #"Salinity_InterstitialWater"#commandArgs()[7]
o <- commandArgs()[8]
## ##
df <- readRDS(f)
dfc <- MRcounts(df,norm=T)
dfc.t <- t(MRcounts(df,norm=T))

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

shplot <- ggplot(index,aes(x=factor(index[[assocvar]]),y=shannon)) + 
  geom_boxplot(aes(fill=factor(index[[assocvar]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - Shannon index",sep=""))+ylab("Alpha diversity (Shannon entropy)")+
  xlab(assocvar)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave('shannon.pdf',plot=shplot,path=o,width=10,height=8)

simplot <- ggplot(index,aes(x=factor(index[[assocvar]]),y=simpson)) + 
  geom_boxplot(aes(fill=factor(index[[assocvar]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - Simpson index",sep=""))+ylab("Alpha diversity (Simpson)")+
  xlab(assocvar)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave('simpson.pdf',plot=simplot,path=o,width=10,height=8)

invsplot <- ggplot(index,aes(x=factor(index[[assocvar]]),y=invsimpson)) + 
  geom_boxplot(aes(fill=factor(index[[assocvar]]))) + geom_jitter(position=position_jitter(width=0.3))+
  ggtitle(paste("Alpha diversity - inverse Simpson index",sep=""))+ylab("Alpha diversity (inverse Simpson)")+
  xlab(assocvar)+scale_fill_discrete(guide=FALSE)+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title = element_text(lineheight=.12, face="bold"))
ggsave('invsimpson.pdf',plot=invsplot,path=o,width=10,height=8)

pdf(paste(o,"renyi_plot.pdf",sep=''),width=10, height=8)
plot(Ren,main="Renyi diversities")

pdf(paste(o,"rankabund_plot.pdf",sep=''),width=10, height=8)
plot(rankabund,main="Ranked abundance distribution models",pch=20)

colourCount =NROW(index)
base <- colorRampPalette(brewer.pal(8, "Set1"))(colourCount)
if (colourCount <= 9) {palette <- base
}else palette <- c(base,colorRampPalette(brewer.pal(8, "Set3"))(colourCount-9))

pdf(paste(o,"rarefaction_curve.pdf",sep=''),width=10, height=8)
rare <- rarecurve(round(dfc.t), step = 100,col=palette,lwd=3,label=F,ylab="OTUs",xlab="Reads sampled",main="Rarefaction curve")
legend(x=5, y=4,rownames(dfc.t), pch=21,col="#777777", pt.bg=palette,pt.cex=2,cex=.8,bty="n",ncol=9,xjust=0,x.intersp=1,text.width=1515)


## end ##
message("\n ** Index files and plots were succesufy created **\n")
