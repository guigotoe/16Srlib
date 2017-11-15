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
#toolbox <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/age_lib/toolbox.R'
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/age_lib/toolbox.R"
toolbox <- "/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","reshape2","optparse","pheatmap","vegan","clusterSim","Rlof","plyr","igraph","orca"))

## Options ##
#p <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib_test/'
p <- '/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib_test/age/'
r <- paste(p,'coocurrences/',sep='')

option_list <- list(
  make_option(c("-i","--data"),type="character",default=paste(p,'dataFcp.rds',sep=''),#NA,#
              help="Path to input rds file"),
  make_option(c("-o","--out"),type="character",default=paste(r,'2017/',sep=''),#
              help="Path to output directory [default %default]"),
  make_option(c("-s","--shared"),type="double",default=100,
              help="OTUs shared by the samples: 0-1; default: %default"),
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
pData(df.r) <- pData(df.r)[order(pData(df.r)$Age),]  # sorting phenotype according Age
df.r <- df.r[,rownames(pData(df.r))]                 # sorting samples according Age
pData(df.r)$LLI <-rep('normal',NROW(pData(df.r)))
pData(df.r)$LLI[pData(df.r)$Age>90] <- "lli"
pData(df.r)$LLI <- as.factor(pData(df.r)$LLI)
pData(df.r)$group <- as.factor(unlist(apply(pData(df.r),1,function(x){
  if (as.numeric(x[["Age"]]) < 40) {return("G1")
  }else if (as.numeric(x[["Age"]]) > 39 & as.numeric(x[["Age"]]) < 60){ return("G2")
  }else if (as.numeric(x[["Age"]]) > 59 & as.numeric(x[["Age"]]) < 80){ return("G3")
  }else if (as.numeric(x[["Age"]]) > 79 & as.numeric(x[["Age"]]) < 90){ return("G4")
  }else if (as.numeric(x[["Age"]]) > 79) return("LLI")
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
phenotype <- list()
#phenotype <- as.data.frame(matrix(data=NA,nrow=0,ncol=length(colnames(pData(df)))))
#colnames(phenotype) <- colnames(pData(df))
j="males"
i <- "LLI"
df <- males
### by age group ###
title <- paste(j,i,sep="_")
samplesToKeepByGroup <- which(pData(df)$group==i)
#df.g <- df.g
x <- 2
## males: G1=1.3 G2=1.3, G3=1.5, G4=1.3; LLI - crush
## females: G1=1.5 G2=1.1 G3=1.4, G4=1.4; LLI=2
df.g <- df[,samplesToKeepByGroup]
## remove outliers using Local Outlier Factor (LOF) approach##
dfc <- t(MRcounts(df.g,norm=T,log=T))
distmat <- vegdist(dfc,"bray")
nmds <- metaMDS(distmat,trace=F)
ordiplot(nmds,type="p",display="sites")
sc <- scores(nmds)
outliers <- which(lof(sc,nrow(dfc)/3)>x)#
plot(lof(sc,nrow(dfc)/3))
abline(h=x,col = "gray60", lty = 3)
if(length(outliers)!=0) df.g <- df.g[,-outliers]
phenotype[[title]] <- pData(df.g)
## end outliers##
## coocurrence for each cluster ##
source(toolbox)
net <- net.mb(df.g,nc=2,tl=df.f)
Nets[[title]] <- net
knet <- net$graph
plot(net$graph, layout=net$layout, vertex.size=net$vsize, vertex.label=NA, main=paste(j,i,sep=" "),edge.width=2)#abs(E(net$graph)$weight)*40)
par(mfrow=c(1,1))
pdf(paste(opt$out,title,"_net.pdf",sep=''),width=8, height=5)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$pcol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$plev, title="Phylum",pch=21,col="#777777", pt.bg=knet$pcol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$ccol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$clev, title="Class",pch=21,col="#777777", pt.bg=knet$ccol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$fcol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$flev, title="Family",pch=21,col="#777777", pt.bg=knet$fcol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$gcol,vertex.label=NA)
legend(x=-2.2, y=-1.1,knet$glev,title="Genus",pch=21,col="#777777", pt.bg=knet$gcol, pt.cex=1.5, cex=.6, bty="n", ncol=5)
#legend(x=1.1, y=1.3,knet$glev,title="Genus",pch=21,col="#777777", pt.bg=knet$gcol, pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

#### all data together ###
saveRDS(Nets,file=paste(opt$out,'Nets.rds',sep=''))
saveRDS(phenotype,file=paste(opt$out,'Clean_phenotype.rds',sep=''))
ph <- ldply(phenotype, data.frame)
ph <- ph[order(ph$Age),]
Nets.a <- list()
phenotype.a <- list()

j="All"
i <- "ALL"
df <- df.f
### by age group ###
title <- paste(j,i,sep="_")
samplesToKeepByGroup <- which(pData(df)$group==i)
df.g <- df
x <- 1.4
## All: G1=1.4 G2=1.2, G3=1.4, G4=1.1; LLI=1.2; All=1.4
df.g <- df[,samplesToKeepByGroup]
## remove outliers using Local Outlier Factor (LOF) approach##
dfc <- t(MRcounts(df.g,norm=T,log=T))
distmat <- vegdist(dfc,"bray")
nmds <- metaMDS(distmat,trace=F)
sc <- scores(nmds)
l <- round(lof(sc,nrow(dfc)/2),1)
outliers <- which(l>x)
cols <- rep("grey",length(l))
cols[outliers]="red"
fig <- ordiplot(nmds,type="none",display="sites")
points(fig,"sites",pc=19,col=cols,cex=0.9)
text(sc,labels=l,pos=3,cex=0.7)
#plot(lof(sc,nrow(dfc)/3))
if(length(outliers)!=0) df.g <- df.g[,-outliers]
phenotype.a[[title]] <- pData(df.g)
## end outliers##
## coocurrence for each cluster ##
source(toolbox)
net <- net.mb(df.g,nc=2,tl=df.f)
Nets.a[[title]] <- net
knet <- net$graph
i <- 'ALL'
ss <- paste('All',i,sep='_')
title <- paste(j,i,sep="_")
net <- Nets.a[[ss]]
knet <- Nets.a[[ss]]$graph
#net$layout <- layout.fruchterman.reingold(net$graph)
plot(net$graph, layout=net$layout, vertex.size=net$vsize, vertex.label=NA, main=paste(j,i,sep=" "),edge.width=2)#abs(E(net$graph)$weight)*40)
pdf(paste(opt$out,title,"_net.pdf",sep=''),width=8, height=5)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$pcol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$plev, title="Phylum",pch=21,col="#777777", pt.bg=knet$pcol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$ccol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$clev, title="Class",pch=21,col="#777777", pt.bg=knet$ccol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$fcol,vertex.label=NA)
legend(x=1.1, y=1.3,knet$flev, title="Family",pch=21,col="#777777", pt.bg=knet$fcol, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$gcol,vertex.label=NA)
legend(x=-2.2, y=-1.1,knet$glev,title="Genus",pch=21,col="#777777", pt.bg=knet$gcol, pt.cex=1.5, cex=.6, bty="n", ncol=5)
#legend(x=1.1, y=1.3,knet$glev,title="Genus",pch=21,col="#777777", pt.bg=knet$gcol, pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()
levels <- taxon.level(df.f,"Family")
partition_assignments <- unlist(lapply(V(knet)$Family, function(x) which(levels==x)))
knet$fassort <- assortativity_nominal(knet,partition_assignments, directed=F)
## Community detection
ncomunities <- length(sizes(knet$ceb)[sizes(knet$ceb)>3])
com <- sizes(knet$ceb)[sizes(knet$ceb)>3]
mark <- lapply(names(com),function(x) which(membership(knet$ceb)==x))
markcol <- adjustcolor(colorRampPalette(brewer.pal(8, "Greys")[4:3])(ncomunities),alpha.f=0.4)
## k-core; The k-core is the maximal subgraph in which every node has degree of at least k.
pdf(paste(opt$out,title,"_net_community.pdf",sep=''),width=8, height=5)
plot(knet$ceb,knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$pcol,vertex.label=NA)
#legend(x=1.1, y=1.3,knet$plev, title="Phylum",pch=21,col="#777777", pt.bg=knet$pcol, pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=-1.2,labels=paste("No. communities > 3 members: ",ncomunities,sep=""),pos=2)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=net$vsize, vertex.color=V(knet)$fcol,vertex.label=NA,
     mark.groups=mark,mark.col=markcol, mark.border="black")
legend(x=1.1, y=1.3,knet$flev, title="Family",pch=21,col="#777777", pt.bg=knet$fcol, pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=-1.2,labels=paste("Homophily: ",round(knet$fassort,2),sep=""),pos=2)
plot(knet, layout=net$layout,rescale=T,edge.curved=.15,edge.width=2,main=paste(j,i,sep=" "),
     vertex.size=knet$kc*5, vertex.color=V(knet)$fcol,vertex.label=NA,
     mark.groups=mark,mark.col=markcol, mark.border="black")
legend(x=1.1, y=1.3,knet$flev, title="Family",pch=21,col="#777777", pt.bg=knet$fcol, pt.cex=2, cex=.8, bty="n", ncol=1)
text(x=-1.2,labels=paste("degree max.: ",max(knet$kc),sep=""),pos=2)
dev.off()

## end coocurrence ##
saveRDS(Nets.a,file=paste(opt$out,'Nets.a.rds',sep=''))
saveRDS(phenotype.a,file=paste(opt$out,'Clean_phenotype.a.rds',sep=''))
Nets.a <- readRDS('/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/age/Nets.a.rds')
net1 <- Nets.a$All_G1$graph
net2 <- Nets.a$All_G2$graph
net3 <- Nets.a$All_G3$graph
net4 <- Nets.a$All_G4$graph
net5 <- Nets.a$All_LLI$graph
net6 <- Nets.a$All_ALL$graph
plot(net1, layout=Nets.a$All_G1$layout, vertex.size=Nets.a$All_G1$vsize, vertex.label=NA,,edge.width=2)#abs(E(net$graph)$weight)*40)

ddnet1 <- degree_distribution(net1,cumulative=F)
ddnet2 <- degree_distribution(net2,cumulative=F)
ddnet3 <- degree_distribution(net3,cumulative=F)
ddnet4 <- degree_distribution(net4,cumulative=F)
ddnet5 <- degree_distribution(net5,cumulative=F)
ddnet6 <- degree_distribution(net6,cumulative=F)

avedeg1 <- sum(seq_along(ddnet1)*ddnet1)-1
avedeg2 <- sum(seq_along(ddnet2)*ddnet2)-1
avedeg3 <- sum(seq_along(ddnet3)*ddnet3)-1
avedeg4 <- sum(seq_along(ddnet4)*ddnet4)-1
avedeg5 <- sum(seq_along(ddnet5)*ddnet5)-1
avedeg6 <- sum(seq_along(ddnet6)*ddnet5)-1

pdf(paste(opt$out,"ALL_DegreeDistribution.pdf",sep=''),width=8, height=5)
plot(seq_along(ddnet1)-1,ddnet1,type='b',xlim=c(0,15),ylim=c(0,0.5),ylab="Frequency",xlab="Degree",col=brewer.pal(8, "Set1")[1])
points(seq_along(ddnet2),ddnet2,type='b',col=brewer.pal(8, "Set1")[2])
points(seq_along(ddnet3),ddnet3,type='b',col=brewer.pal(8, "Set1")[3])
points(seq_along(ddnet4),ddnet4,type='b',col=brewer.pal(8, "Set1")[4])
points(seq_along(ddnet5),ddnet5,type='b',col=brewer.pal(8, "Set1")[5])
points(seq_along(ddnet6),ddnet6,type='b',col=brewer.pal(8, "Set1")[7])
legend("topright",c("G1","G2","G3","G4","LLI  ","ALL"),col=brewer.pal(8, "Set1")[c(1:5,7)],pch=1,lty=1)
dev.off()

natcon <- function(ig) {
  N   <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc  <- log(mean(exp(evals)))
  nc / (N - log(N))
}
ig <- net3
nc.attack <- function(ig) {
  E(ig)$weight <- abs(E(ig)$weight)
  hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*0.8), function(i) {
    ind <- hubord[1:i]
    #print(length(hubord))
    tmp <- delete_vertices(ig, V(ig)[[ind]])
    natcon(tmp)
  }) }
nc.g1 <- nc.attack(net1)
nc.g2 <- nc.attack(net2)
nc.g3 <- nc.attack(net3)
nc.g4 <- nc.attack(net4)
nc.g5 <- nc.attack(net5)
nc.g6 <- nc.attack(net6)

pdf(paste(opt$out,"ALL_NetworkStability.pdf",sep=''),width=8, height=5)
plot(seq(0,.8,len=length(nc.g1)), nc.g1, type='l', ylim=c(0,0.013),xlab="Proportion of removed nodes", ylab="natural connectivity")
points(seq(0,.8,len=length(nc.g2)), nc.g2, type='l',col=brewer.pal(8, "Set1")[2])
points(seq(0,.8,len=length(nc.g3)), nc.g3, type='l',ylim=c(0,max(nc.g1)),col=brewer.pal(8, "Set1")[3])
points(seq(0,.8,len=length(nc.g4)), nc.g4, type='l',col=brewer.pal(8, "Set1")[4])
points(seq(0,.8,len=length(nc.g5)), nc.g5, type='l',col=brewer.pal(8, "Set1")[5])
points(seq(0,.8,len=length(nc.g6)), nc.g6, type='l',col=brewer.pal(8, "Set1")[7])
legend("topright",c("G1","G2","G3","G4","LLI  ","ALL"),col=brewer.pal(8, "Set1")[c(1:5,7)],pch=1,lty=1)
dev.off()

hist(diff(nc.g1), breaks=30, col='red', ylim=c(0,80), main="Fragility rates")


selist <- list(Nets.a$All_G1$mb, Nets.a$All_G2$mb,Nets.a$All_G3$mb, Nets.a$All_G4$mb,Nets.a$All_LLI$mb, Nets.a$All_ALL$mb)
orbits <- c(0,2,5,7,8,10,11,6,9,4,1) + 1
graphlets <- lapply(selist, function(x) count4(summary(Matrix::triu(x$refit))[,1:2])[,orbits])

###
## making a movie ##
packages(c("network","sna","ergm","tsna","ndtv","visNetwork","networkDynamic"))
netlist <- list()
for (i in names(Nets.a)){
  print(i)
  netlist[[i]] <- as_adj(Nets.a[[as.character(i)]]$graph,edge=T,sparse=T)
}
netlist<-lapply(netlist,as.network.matrix,matrix.type='adjacency')
names(netlist)
allnets<-networkDynamic(network.list=netlist)
p.names <- V(Nets.a[[as.character(i)]]$graph)$Phylum
c.names <- V(Nets.a[[as.character(i)]]$graph)$Class
f.names <- V(Nets.a[[as.character(i)]]$graph)$Family
g.names <- V(Nets.a[[as.character(i)]]$graph)$Genus
network.vertex.names(allnets)<-c.names
#mnet%v%'family' <- f.names
#mnet%v%'phyllum' <- p.names

allnets%n%'net.obs.period'<-list(
  observations=list(c(0,length(names(allnets)))),
  mode="discrete", 
  time.increment=1,
  time.unit="book volume")

render.animation(allnets)
compute.animation(allnets,
                  animation.mode='MDSJ',
                  default.dist=2,
                  verbose=FALSE)
slice.par <- list(start = 1, end = 5, interval = 1,
                  aggregate.dur = 1, rule = "latest")
compute.animation(allnets,
                  animation.mode='MDSJ',
                  default.dist=2,
                  verbose=FALSE)
list.vertex.attributes(allnets)
render.animation(allnets)

render.d3movie(nt,
               render.par=list(tween.frames=30,show.time=F),
               vertex.cex=0.8,label.cex=0.8,label.col='gray',
               displaylabels=F,
               # make shape relate to school year
               #vertex.sides=mnet%v%'schoolyear'-1983,
               # color by gender
               #vertex.col=ifelse(harry_potter_support%v%'gender'==1,'blue','green'),
               edge.col="darkgray",
               vertex.tooltip=function(slice){paste('name:',slice%v%'vertex.names','<br>',
                                                    'status:', slice%v%'testatus')},
               output.mode = 'htmlWidget'
)
m <- mnet
par(mfrow=c(1,1))
par(mar = c(3, 3, 2, 4))

proximity.timeline(allnets,default.dist=6,
                   mode='sammon',labels.at=124,vertex.cex=4)
plot(tEdgeFormation(allnets))
plot(tSnaStats(allnets,'gtrans'))
plot(tErgmStats(allnets,'meandeg'),main='Mean degree of Men nets')

timeline(allnets,slice.par=list(start=0,end=25,interval=1,aggregate.dur=1,rule='latest'),
         plot.vertex.spells=FALSE)

timePrism(allnets,at=c(1,10,20),displaylabels=T,planes=T,label.cex=0.5)
timeline(allnets)



###


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

