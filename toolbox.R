####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: September 2016
# Created: 29 October 2015
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
#
# toolbox.R this script contains several functions
# commonly used by the 16Srlib pipeline for 16S analysis.
#
####################################################

packages <- function(requirements){
  has   <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    #setRepositories(ind=1:10)
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has],repos=getOption("repos")[!is.na(getOption("repos"))]) #repos="https://cran.uni-muenster.de/"
  }
  suppressPackageStartupMessages(lapply(requirements,require,character.only=T,quietly=T))
}

numDFtranspose <- function(df){
  x <- t(df)
  colnames(x) <- x[1,]
  x <- x[-1,]
  rowid <- rownames(x)
  z <- as.data.frame(apply(x, 2, as.numeric))
  rownames(z) <- rowid
  return(z)
}

changename <- function(df,col){
  df[[col]] <- gsub("unclassified",'',df[[col]])
  df[[col]] <- gsub("",'',df[[col]])
  return(df)
}

mothur.taxonomy <- function(taxonomy){
  packages(c("reshape2"))
  taxnames <- c('Kingdom','Phylum','Class','Order','Family','Genus')
  bootsraps <- c('bootstrap1','bootstrap2','bootstrap3','bootstrap4','bootstrap5','bootstrap6')
  x <- colsplit(taxonomy$Taxonomy,pattern="\\([[:digit:]]*\\);",names = taxnames)
  x <- apply(x,2,function(x) as.factor(gsub("\\([[:digit:]]*\\);",'',x)))
  y <- colsplit(taxonomy$Taxonomy,pattern="\\;[[[:alnum:]]_-]*",names = bootsraps)
  y <- apply(y,2,function(y) as.double(gsub("[Bacteria \\(\\)\\;]",'',y)))
  #for (i in taxnames){
  #  x <- changename(x,i)
  #}
  #x[x==""|is.na(x)]  <- 'unclassified'
  df <- cbind(taxonomy[,1:2],x,y)
  rownames(df) <- df[,1]
  return(df)
}

mothur.counts <- function(counts){
  counts$X <- NULL
  counts <- counts[,-c(1,3)]
  colnames(counts)[1] <- "id"
  df <- numDFtranspose(counts)
  return(df)
}

mothur.metadata <- function(metadata){
  if (is.na(table(duplicated(metadata[,1]))["TRUE"])){row.names(metadata) <- metadata[,1]
  } else {rownames(df) <- make.names(metadata[,1], unique=TRUE)}
  return(metadata)
}

#biomf <- opt$biom
#metadataf <- opt$metadata
mothur.biom <- function(biomf,metadataf){
  ## biom = path, metadata = path
  packages(c("biomformat","metagenomeSeq"))
  biom_file <- read_biom(biomf)
  bdata <- biom2MRexperiment(biom_file)
  taxnames <- c('Kingdom','Phylum','Class','Order','Family','Genus','Specie')
  colnames(fData(bdata)) <- taxnames
  mdata <- mothur.metadata(read.table(metadataf,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
  ord = match(colnames(MRcounts(bdata)),rownames(mdata))
  mdata = mdata[ord,]
  if(length(colnames(MRcounts(bdata)))!=length(union(colnames(MRcounts(bdata)),rownames(mdata)))){
    message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
    quit()
  }
  taxdf <- fData(bdata)
  taxdf$Size <- apply(MRcounts(bdata),1,sum)
  counts <- as.data.frame(MRcounts(bdata))
  counts$Size <- as.numeric(unlist(lapply(rownames(MRcounts(bdata)),function(x) taxdf[x,'Size'])))
  ##* new names.. better not because picrust use the old ones *##
  #taxdf <- taxdf[with(taxdf,order(-Size)),]
  #counts <- counts[with(counts,order(-Size)),]
  #taxdf$OTU <- paste('GG_OTU',1:nrow(taxdf),sep='')
  #rownames(counts) <- taxdf$OTU
  #rownames(taxdf) <- taxdf$OTU
  data.mr <- newMRexperiment(counts[,1:(ncol(counts)-1)],phenoData=AnnotatedDataFrame(mdata),featureData=AnnotatedDataFrame(taxdf[,c(8,1:7)]))#
  #data.mr <- newMRexperiment(MRcounts(bdata),phenoData=AnnotatedDataFrame(mdata),featureData=AnnotatedDataFrame(fData(bdata)))
  return(data.mr)
}

mothur.usingcounts <- function(countsf,metadataf,taxonomyf){
  packages('metagenomeSeq')
  metadata <- mothur.metadata(read.table(metadataf,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
  taxonomy <- mothur.taxonomy(read.table(taxonomyf,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
  counts <- mothur.counts(read.table(countsf,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
  #taxonomy <- subset(taxonomy,taxonomy$Kingdom!="unclassified"& taxonomy$Phylum!="unclassified")
  #counts["7751662AGE_a_G4"] <- unlist(apply(counts[c("7344072AGE1_a_G4","6897722AGE_a_G4","662586SPC_a_G4","7344072AGE1_a_G4","9457743AGE1_a_G4")][,,],1,mean))
  ord = match(colnames(counts),rownames(metadata))
  metadata = metadata[ord,]
  if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
    message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
    quit()
  }
  metadata$libsize <- apply(counts,2,sum)
  data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
  return(data)
}

libsizefilter <- function(data,libt,outpath){
  outliers <- boxplot.stats(pData(data)$libsize)$out
  bottomth <- mean(pData(data)$libsize[-which(pData(data)$libsize%in%outliers)],na.rm=T)-(as.numeric(libt)*sd(pData(data)$libsize[-which(pData(data)$libsize%in%outliers)],na.rm=T))
  pdf(paste(outpath,'Library size',sep=''), width=10, height=5)
  plot(pData(data)$libsize,col=ifelse(pData(data)$libsize>bottomth,"blue","red"),pch=20,
     xlab="Sample", ylab="No. of reads")
  abline(h=bottomth,col="RED",lty=2)
  dev.off()
  keepSamples <- pData(data)[pData(data)$libsize>bottomth,1]
  return(data[,keepSamples])
}

pipts.otu <- function(p.counts,p.metadata){
  metadata <- mothur.metadata(read.table(p.metadata,header=T,sep="\t",na.strings=c("","NA")))
  df <- read.table(p.counts,header=T,sep="\t",quote='',row.names=1,as.is=T)
  colnames(df) <- gsub("^X", "",  colnames(df))
  counts <- df[,-which(colnames(df)=="taxonomy")]
  taxnames <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  x <- colsplit(df$taxonomy,pattern="; ",names = taxnames)
  row.names(x) <- row.names(counts)
  bootsraps <- c('bootstrap1','bootstrap2','bootstrap3','bootstrap4','bootstrap5','bootstrap6','bootstrap7')
  y <- data.frame(matrix(NA, nrow=length(row.names(x)), ncol=length(bootsraps),dimnames=list(row.names(x),bootsraps)))
  taxonomy <- cbind(x,y)
  ord = match(colnames(counts),rownames(metadata))
  metadata = metadata[ord,]
  if(length(colnames(counts))!=length(union(colnames(counts),rownames(metadata)))){
    message("**Error: Count-sample's names don't match with Metadata-sample's names!\n\t*Check the files and try again\n")
    quit()
  }
  taxonomy$OTU <- rownames(taxonomy)
  taxonomy$Size <- apply(counts,1,sum)
  taxonomy <- taxonomy[,c(15,16,1:14)]
  metadata$libsize <- apply(counts,2,sum)
  data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
  return(data)
}
  #opath <- opt$out
replicas.analysis <- function(data,opath){
  packages(c('metagenomeSeq','gridExtra'))
  metadata <- pData(data)
  techrep <- metadata$CC[duplicated(metadata$CC)]
  message('normalizing data...')
  counts.nl <- MRcounts(data,norm=T,log=T)
  xyi=c()
  plots <- list()
  message('correlating replicates...')
  for (i in seq(length(techrep))){
    vec <- 1:nrow(metadata[metadata$CC==techrep[i],])
    combs <- combn(vec,2)
    for (j in 1:ncol(combs)){
      df.s <- metadata[metadata$CC==techrep[i],][combs[,j],]
      y <- counts.nl[,df.s$ID[1]]
      x <- counts.nl[,df.s$ID[2]]
      df <- data.frame(x=x,y=y)
      df$ratio <- unlist(apply(df,1,function(x){
        if(is.finite(x[1]/x[2])){if(x[1]/x[2]==0)return(x[2]) else return(x[1]/x[2])
        }else if(is.nan(x[1]/x[2])){return(1)
        }else if(is.infinite(x[1]/x[2])) return(x[1])
      }))
      #plot(density(df$ratio))
      corrxy <- cor.test(df$x,df$y)
      if(corrxy$estimate>0.8){
        minval <- mean(df$ratio)-sd(df$ratio)
        maxval <- mean(df$ratio)+sd(df$ratio)
        df$highlight <- unlist(lapply(df$ratio,function(x) if(x!=1){return("highlight")}else{return("normal")}))
        dz <- subset(df,df$highlight=="highlight")
        z <- abs(dz$x-dz$y) # how far is one value to another.
        xyintercept <- quantile(z,probs=0.95)
        #hist(z)
        #abline(v=xyintercept)
        #xyintercept <- max(z)
        xyi <- c(xyi,xyintercept)
        df$keept <- unlist(apply(df,1,function(x) if(as.numeric(x[1])>xyintercept&as.numeric(x[2])>xyintercept){return("keept")}else{return("drop")}))
        print(c(paste(i,j,sep="_"),corrxy$estimate,xyintercept,'OTUs_retained'=length(df$keept[df$keept=='keept'])))
        mycolours <- c("keept" = "red", "drop" = "black")
        dz2 <- subset(df,df$highlight=="highlight")
        fig <- ggplot(dz2,aes(x=x,y=y))+geom_point(aes(alpha=1/100,colour=keept))+
          geom_hline(yintercept=xyintercept,co)+geom_vline(xintercept=xyintercept)+scale_colour_manual(values=mycolours)+
          ggtitle(paste(i,'_',j," Pcor: ",round(corrxy$estimate,2),' Int.val(log2): ',round(xyintercept,2),' Int.val(counts): ',round(2^xyintercept,2),sep=""))+
          theme(legend.position="none",plot.title=element_text(lineheight=.8, face="bold"))
        plots[[paste(i,j,sep="_")]] <- fig
        #ggsave(filename=paste(opt$out,i,'_replicates.pdf',sep=''),device="pdf")
      }
    }
    plotfiles <- c()
    pdf(paste(opath,'TechRep_scatter.pdf',sep=''),onefile=T,width=9,height=7)
    for (q in seq(length(plots))){do.call("grid.arrange",plots[q])}
    dev.off()
  }
  #plotDensities(xyi)
  ##threshold <- quantile(xyi,probs=0.95)
  #shapiro.test(xyi)
  #hist(xyi)
  #abline(v=threshold)
  #abline(v=(mean(xyi)+sd(xyi)))0
  #round(2^mean(threshold))
  #summary(xyi)
  #round(2^threshold)
  #round(2^mean(xyi))
  #round(2^(mean(xyi)+sd(xyi)))
  otu.coverage <- round(2^(mean(xyi)+sd(xyi))) # threshold implemented as mean+sd. Stringent enough infered from data.
  ###otu.coverage <- round(2^mean(5.5)) ## manual adjust
  ##otu.coverage <- 6 ##*** because I did prev. manually and is log=T then i should use xyi value: manualy placed as 6
  #otu.coverage <- round(2^mean(xyi))
  ## removing replications
  rep <- metadata[metadata$CC%in%techrep,]
  samplesToKeep <- rownames(metadata)[!rownames(metadata)%in%rownames(rep[rep$rep=="a",])]
  data.f1 <- data[,samplesToKeep]
  remove(counts.nl)
  return(list(data=data.f1,coverage=otu.coverage,interceptcoords=xyi))
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

dismatplot <- function(sampleDists,samplenames,title,o){
  packages(c("pheatmap"))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- samplenames #paste(rownames(pData(df)), sep="" )
  colnames(sampleDistMatrix) <- samplenames
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
  pdf(paste(o,"Heatmap_",gsub(" ",'_',title),".pdf",sep=""),pointsize=12, width=10,height=7,onefile=FALSE)
  pheatmap(sampleDistMatrix,main=title,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
}

taxon.level <- function(df,tlevel,pmin=63){
  packages("matrixStats")
  df <- aggTax(df,lvl=tlevel) 
  ma <- rowMeans(df) # Mean abundance of this taxa
  p <- apply(MRcounts(df),1,function (x) sum(x>63)/dim(df)[2]) # Total presence = presence/No. Individuals => max = 1
  taxorder <- ma*p
  taxorder <- sort(taxorder,decreasing=T)
  return(names(taxorder))
}
#tl <- 'Phylum'

taxprop <- function(dfp,v,tl,o,limit=15,u=T,z=T,stack=F,legendncol=NA,textsize=NA){
  ## input: 1) MRexperiment (proportions,features,metadata) 2) reference variable: x axis 3) taxlevel 4) output path
  ## u=F -> with unclisified -> z=F -> with othersas.numeric(df.t[i,])>=quantile(as.numeric(df.t[i,]), probs=0.9)
  ## limit 0.4 -> selection of Taxa which comprises 40% of the proportion at each sample.
  ## sum columns with same name -- OTUs with same taxonomical level:
  #dfxta <- as.data.frame(sapply(unique(names(dfxt)[duplicated(names(dfxt))]), function(x) Reduce("+", dfxt[ , grep(x, names(dfxt))])))
  packages(c("matrixStats","RAM"))
  #tl="Phylum"
  title <- paste(tl," abundance distribution",sep='')
  if (!"ID"%in%colnames(pData(dfp)))pData(dfp)$ID <- pData(dfp)[,1]
  df <- aggTax(dfp,lvl=tl)
  df.t<- as.data.frame(t(MRcounts(df)))
  ##########################################
  ## * Selecting x (limit) OTUs to plot *##
  OTUsToKeep <- c()
  byfile_OTUabundances <- data.frame(OTU=colnames(df.t)) ## To export later all the abundances according files
  for (i in 1:nrow(df.t)){ 
    ## Keeps OTUs according limit -- so if limit is >1 represent the maximum number of OTUs to plot (most abundant)
    ## if limit <1 then will keep OTUs based on the abundance proportion. By default limit = 15 
    pa <- data.frame(OTU=colnames(df.t[i,]),freq=as.numeric(df.t[i,]))
    pa <- pa[with(pa,order(-freq)),]
    pa$Fcum <- cumsum(pa$freq)
    rownames(pa) <- 1:nrow(pa)
    byfile_OTUabundances[[rownames(df.t)[i]]] <- pa[match(byfile_OTUabundances$OTU, pa$OTU),]$freq
    if(!u){
      if(any(as.character(pa$OTU) == 'Bacteria_unclassified')) pa <- pa[-which(pa$OTU=='Bacteria_unclassified'),]
    }
    if (limit > 1){palim <- pa[1:limit,]
    }else {palim <- pa[pa$Fcum<=limit,]; if(nrow(palim)<=10) palim <- pa[1:10,]} # for phylum mainly but if doesnt work for others
    OTUsToKeep <- c(OTUsToKeep,as.character(palim$OTU))
  }
  OTUsToKeep <- unique(OTUsToKeep)
  ## However the prev. step we will end with more OTUs that expected. We should stick to our LIMIT then:
  # We will use the top x (LIMIT) according to the average abundance across samples. 
  byfile_OTUabundances$rowMean <- rowMeans(byfile_OTUabundances[2:ncol(byfile_OTUabundances)])
  OTUsToKeep.filter <- byfile_OTUabundances[byfile_OTUabundances$OTU%in%OTUsToKeep,c(1,ncol(byfile_OTUabundances))]
  OTUsToKeep.filter <- OTUsToKeep.filter[with(OTUsToKeep.filter,order(-rowMean)),]
  OTUsToKeep <- OTUsToKeep.filter[1:limit,]$OTU
  #####################################################
  
  df.tk <- df.t[,which(colnames(df.t)%in%OTUsToKeep)] ## keeping taxon with high proportion according 'limit'
  df.tk$Others <- unlist(apply(df.tk,1,function(x){1-sum(x)}))
  
  ###** Writing df.tk and byfile_OTUabundances:
  write.table(round(df.tk,2),file=paste(o,gsub(" ",'_',title),'_',v,'_',ncol(df.tk),'OTU.txt',sep=''),sep="\t",quote=F,row.names=F)
  byfile_OTUabundances_x <- round(byfile_OTUabundances[2:ncol(byfile_OTUabundances)],5)
  byfile_OTUabundances_x$OTU <- byfile_OTUabundances$OTU
  byfile_OTUabundances_x <- byfile_OTUabundances_x[with(byfile_OTUabundances_x,order(-rowMean)),]
  byfile_OTUabundances_x$rMeanCum <- cumsum(byfile_OTUabundances_x$rowMean)
  write.table(byfile_OTUabundances_x,file=paste(o,gsub(" ",'_',title),'_',v,'.txt',sep=''),sep="\t",quote=F,row.names=F)
  ###
  
  gm <- cbind(df.tk,pData(dfp)[v],ID=pData(dfp)$ID) 
  gm_m <- melt(gm,id.vars=c(v,"ID"))
  colnames(gm_m)[which(colnames(gm_m)=='variable')] <- 'taxa'
  ##################################################################
  #*** following line is a depreciated process *** not used anymore
  #Because some taxa now are xOthers so we need to summ their values
  #gmx2 <- aggregate(as.formula(paste("value~",'ID+taxa+',paste(v,collapse="+"))),data=gm_m,FUN=sum) 
  ##################################################################
  gmx2 <- gm_m
  taxorder <- data.frame(taxa=levels(gmx2$taxa),rate=rep(0,length(levels(gmx2$taxa))))
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
  
  ut <- ''
  ot <- ''
  if(!u){gmx2 <- subset(gmx2,gmx2$taxa!='unclassified');ut <- '_Xuncl'}
  if(!z){gmx2 <- subset(gmx2,gmx2$taxa!='Others');ot <- '_Xothers'}
  
  #pointsplot <- 
  
  #geom_point(aes(colour=gmx2$taxa,fill=gmx2$taxa),position=position_jitter(width=0.22),alpha = 0.9)+
  #violinplot <- ggplot(gmx2,aes(x=as.factor(gmx2[[v]]),y=value,fill=taxa))+  
  #  geom_violin(aes(colour=gmx2$taxa,fill=gmx2$taxa))+
  #  geom_point(aes(colour=gmx2$taxa,fill=gmx2$taxa,y=mean(gmx2$value)),color="black", size = 2) + 
  #  geom_errorbar(aes(colour=gmx2$taxa,y=mean(gmx2$value),ymin=mean(gmx2$value)-sd(gmx2$value),ymax=mean(gmx2$value)+sd(gmx2$value)), 
  #                color = "black", width = 0.2)+ 
  #  scale_fill_manual(name=tl,values=palette)+
  #  scale_colour_manual(name=tl,values=palette)+
  #  scale_x_discrete(v)+ylab("Proportion")+
  #  ggtitle(paste(title,sep=""))+
  #  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12,face="bold"),
  #        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
  #        legend.position="bottom",legend.box="horizontal",
  #        legend.text = element_text(size=9),
  #        legend.title = element_text(size=10, face="bold"),
  #        plot.title = element_text(lineheight=.12, face="bold"))
  if (is.na(textsize)){
    if(length(levels(gmx2$taxa))<=14){textsize=10;keysize=1}
    if(length(levels(gmx2$taxa))>14){textsize=8;keysize=0.8}
    #if(length(levels(gmx2$taxa))>26) {textsize=6.5;keysize=0.7}
    if(length(levels(gmx2$taxa))>27) {textsize=6;keysize=0.6}
  }
  
  distplot <- ggplot(gmx2,aes(x=as.factor(gmx2[[v]]),y=value,fill=taxa))+
    geom_boxplot(lwd=0.2,width=0.75)+#geom_jitter(position=position_jitter(width=.2), size=0.5)+
    scale_fill_manual(name=tl,values=palette) +
    scale_x_discrete(v,expand=c(0,0.4))+ylab("Proportion")+
    ggtitle(paste(title,sep=""))+
    #guides(shape=guide_legend(override.aes = list(size=20)))+a
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=20,face="bold"),
          axis.text.y  = element_text(size=20),
          axis.title.y=element_text(size=20,face="bold"),
          axis.title.x=element_blank(),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          legend.position="bottom",legend.box="horizontal",
          legend.text = element_text(size=textsize),
          legend.title = element_text(size=10, face="bold"),
          plot.title = element_text(lineheight=.12, face="bold"),
          legend.key.size=unit(keysize,'line'),
          plot.background= element_blank())
  if(!is.na(legendncol))distplot <- distplot+guides(fill=guide_legend(ncol=legendncol))
  #distplot
  stackplot <- ggplot(gmx2,aes(x=as.factor(gmx2[[v]]),y=value,fill=taxa))+geom_bar(stat='identity',aes(fill=taxa))+
    scale_fill_manual(name=tl,values=palette) +
    scale_x_discrete(v)+ylab("Proportion")+
    ggtitle(paste(title,sep=""))+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12,face="bold"),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          legend.position="bottom",legend.box="horizontal",
          legend.text = element_text(size=9),
          legend.title = element_text(size=10, face="bold"),
          plot.title = element_text(lineheight=.12, face="bold"))
  #stackplot
  if(stack==T) ploting <- stackplot else ploting <- distplot
  ggsave(paste(gsub(" ",'_',title),ot,ut,'_',v,'.pdf',sep=''),plot=ploting,path=o,width=8,height=6,device="pdf")
  return(byfile_OTUabundances_x)
}

taxonprop <- function(dfp,v,tl,otu,o){
  #not finish yet
  ## input: 1) MRexperiment (proportions,features,metadata) 2) reference variable: x axis 3) taxlevel 4) Taxon name 5) output path
  OTUsToKeep <-  as.character(fData(dfp)[fData(dfp)[[tl]]%in%otu,1])
  dfx <- dfp[OTUsToKeep,]
  df <- aggTax(dfx,lvl=tl)
  
  df.t<- as.data.frame(t(MRcounts(df)))
  gm <- cbind(df.t,pData(df)[v],id=pData(df)[1])
  gm_m <- melt(gm,id.vars=c(v,"ID"))
  gm_m$"taxa" <- gm_m$"variable"
  gmx2 <- aggregate(as.formula(paste("value~",'ID+taxa+',paste(v,collapse="+"))),data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
  title <- paste(tl," abundance distribution",sep='')
  distplot <- ggplot(gmx2,aes(x=as.factor(gmx2[[v]]),y=value,fill=taxa))+geom_boxplot()+#geom_jitter(position=position_jitter(width=.2), size=0.5)+
    #scale_fill_manual(name=tl)# +
    scale_x_discrete(v)+ylab("Proportion")+xlab("taxon")+
    ggtitle(paste(title,sep=""))+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14,face="bold"),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          legend.position="bottom",legend.box="horizontal",
          legend.text = element_text(size=9),
          legend.title = element_text(size=14, face="bold"),
          plot.title = element_text(lineheight=.16, face="bold"))
  ggsave(paste(gsub(" ",'_',title),'_',v,'.pdf',sep=''),plot=distplot,path=o,width=8,height=5,device="pdf")
}

#DEdata <- MRcounts(dfc,norm=T)
#ldata <- MRcounts(dfc,norm=T,log=T)
#design=design
#method=opt$clmethod
#path=opt$out
#prefix=paste(opt$level,'Sign_conf_',paste(opt$conf,collapse='_'),'_',sep="")
#val=opt$clval

getClusters <- function(DEdata,ldata,method=c("PAM","P","K","Km"),design,path,prefix="",w=2,val=NULL,JSD=F,
                        cellwidth=NA,cellheight=NA,text_size=NA,height=8,width=7){
  packages(c("clusterSim","cluster","pheatmap","vegan"))
  if (missing(method))method <- "PAM"
  if (method=="P" & missing(val))return("argument \"val\" is missed")
  # Recomended options: method PAM, val=NULL
  myheatcol = colorpanel(75, 'blue','grey','red')
  data = t(scale(t(as.matrix(ldata)), scale=F)) # center rows by mean substraction
  if (JSD==T){
    gene_dist <- dist.JSD(t(DEdata))
    sample_dist = dist.JSD(DEdata)
    dis <- 'JDS'
  }else{
    gene_dist = dist(data, method='euclidean')
    sample_dist = dist(t(data), method='euclidean')
    dis <- 'Euc'
  }
  hc_genes = hclust(gene_dist, method='complete')
  sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
  sample_dist = dist(t(data), method='euclidean')
  hc_samples = hclust(sample_dist, method='complete')
  
  if (is.null(val)){
    nclusters=NULL
    
    for (k in 1:30) { 
      if (k==1) {
        nclusters[k]=NA 
      } else {
        if (method == "Km"){data.cluster_temp = kmeans(data, centers=k, iter.max=100, nstart=5)$cluster
        }else if (method == "K"){data.cluster_temp = cutree(as.hclust(hc_genes), k=k)
        }else if (method == "PAM"){data.cluster_temp = pam(as.dist(gene_dist),k, diss=TRUE)$clustering
        }
        nclusters[k]=index.G1(data,data.cluster_temp,d = gene_dist,centrotypes="medoids")
      }
    }
    pdf(paste(path,"OptClusterNum_CHindex.pdf",sep=''),width=8, height=5)
    plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
    dev.off()
    val <- which(nclusters==max(nclusters[-1]))[1]
  }
  
  if (method == "Km"){
    # Use K-means clustering to define K gene sets. (use the -K parameter). 
    # This does not leverage the already hierarchically clustered genes as shown in the heatmap, 
    # and instead uses a least-sum-of-squares method to define exactly k gene clusters.
    kmeans_clustering <- kmeans(data, centers=val, iter.max=100, nstart=5)
    gene_partition_assignments = kmeans_clustering$cluster
    file_name <- paste("clusters_fixed_",dis,"_Kmeans_",val,sep="")
  }else if (method == "K"){
    # Cut the hierarchically clustered genes (as shown in the heatmap) into exactly K clusters
    gene_partition_assignments <- cutree(as.hclust(hc_genes), k=val) # cut the tree into k clusters
    file_name <- paste("clusters_fixed_",dis,"_Ktree_",val,sep="")
  }else if (method == "P"){
    #(Recommended) cut the hierarchically clustered gene tree at --Ptree percent height of the tree.
    gene_partition_assignments <- cutree(as.hclust(hc_genes), h=val/100*max(hc_genes$height))
    file_name <- paste("clusters_fixed_",dis,"_P_",val,sep="")
  }else if (method == "PAM"| is.null(method)){
    # Used by Arumugam et al 2011.Partitioning around medoids (PAM) clustering algorithm
    # to cluster the abundance profiles. PAM derives from the basic k-means algorithm, 
    # but has the advantage that it supports any arbitrary distance measure and is more robust 
    # than k-means
    gene_partition_assignments = pam(as.dist(gene_dist),val, diss=TRUE)$clustering
    file_name <- paste("clusters_fixed_",dis,"_PAM_",val,sep="")
  }
  
  if(is.na(text_size)){
    ifelse(nrow(ldata)>50 & nrow(ldata)<=80,yes=text_size <- 5, no=text_size <- 9 )
    if(nrow(ldata)<=10) text_size <- 12
    if(nrow(ldata)>10 & nrow(ldata)<=30) text_size <- 9
    if(nrow(ldata)>30 & nrow(ldata)<=50) text_size <- 7.5
    if(nrow(ldata)>50 & nrow(ldata)<=80) text_size <- 5.5
    if(nrow(ldata)>80 & nrow(ldata)<=110) text_size <- 4.5
    if(nrow(ldata)>110& nrow(ldata)<=130) text_size <- 3.5
    if(nrow(ldata)>130) text_size <-1.3
  }
  
  
  max_cluster_count = max(gene_partition_assignments)
  outdir <- paste(path,prefix,file_name,sep='')
  if(dir.exists(outdir)) message('Out-folder already exist, files will be overwritten')
  dir.create(outdir,showWarnings=F)
  
  ##* Colors *##
  ## set color representation for specific values of the data distribution
  quantile_range <- quantile(ldata, probs = seq(0, 1, 0.2))
  #color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range)-1)
  color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range)+1)
  # colors clusters #
  colourCount = length(unique(gene_partition_assignments))
  #partition_colors = rainbow(colourCount, start=0.4, end=0.95)
  if (colourCount>9){ l <- colourCount }else l <- 9
  partition_colors <- colorRampPalette(brewer.pal(8, "Greys"))(colourCount+1)[-1]#(l)[c(1:colourCount)]
  gene_colors_dframe = data.frame(clusters=gene_partition_assignments, colors=partition_colors[gene_partition_assignments])
  gene_colors_dframe$id <- factor(rownames(gene_colors_dframe),levels=hc_genes$labels[hc_genes$order])
  gene_colors_dframe <- gene_colors_dframe[order(gene_colors_dframe$id),]
  rowann <- data.frame(Cluster=as.factor(gene_colors_dframe[,1]),row.names=rownames(gene_colors_dframe))
  #write.table(gene_colors_dframe, file=paste(path,prefix,file_name,".heatmap.gene_cluster_colors.dat",sep=''), quote=F, sep='  ')  gene_colors = as.matrix(partition_colors[gene_partition_assignments])
  # colors variable#
  varcolours <- length(levels(design[,2])) ## categorical var
  if (varcolours>9){ l <- varcolours }else l <- 9
  sampleColors = colorRampPalette(brewer.pal(8, "Set1"))(l)[c(1:varcolours)]
  sample_colors_dframe = data.frame(samples=levels(design[,2]), colors=sampleColors)
  sample_colors = unlist(lapply(design[,2], function(x) return(as.character(sample_colors_dframe[x,2]))))
  ## anotation list colors for pheatmap
  #ann_colors = list(Cluster=setNames(partition_colors,levels(as.factor(gene_colors_dframe[,1]))))
  ann_colors <- list()
  ann_colors[['Cluster']] <- setNames(partition_colors,levels(as.factor(gene_colors_dframe[,1])))
  ann_colors[[names(design)[2]]] <- setNames(as.character(sample_colors_dframe$colors),sample_colors_dframe$samples)
  ##* end Colors *#
  pdf(paste(path,prefix,file_name,".heatmap.pdf",sep=''),onefile=F,useDingbats=F,height=height,width=width,bg='transparent')
  pheatmap(data[,rownames(design)],col=color_palette,drop_levels=T,fontsize_row=text_size,annotation_names_row=F,
           cluster_rows=hc_genes,annotation_row=rowann,annotation_colors=ann_colors,
           cluster_cols=F,annotation_col=design[,2,drop=F],scale="row",annotation_names_col=F,
           treeheight_col=F,border_color='gray20',cellwidth=cellwidth,cellheight=cellheight)
  dev.off()

  gene_names = rownames(data)
  num_cols = length(data[1,])
  raw_outdir <- paste(outdir,"/rawCounts",sep='')
  dir.create(raw_outdir)
  for (i in 1:max_cluster_count) {
    partition_i = (gene_partition_assignments == i)
    partition_data = data[partition_i,,drop=F]
    #partition_data = partition_data[,hc_samples$order, drop=F]  # Uncomment if you want order samples by hcluster similarities
    partition_data = partition_data[,order(design[,1]),drop=F]
    outfile = paste(outdir,"/",prefix,"cluster_",formatC(i, width=w, flag="0"), "_log_medianCentered.matrix.txt", sep='')
    write.table(partition_data, file=outfile, quote=F, sep="\t")
    raw_p_data = DEdata[partition_i,,drop=F]
    
    #raw_p_data = raw_p_data[,raw_hc_samples$order, drop=F]      # Uncomment if you want order samples by hcluster similarities
    raw_p_data = raw_p_data[,order(design[,1]),drop=F]
    raw_outfile = paste(raw_outdir,"/",prefix,"cluster_",formatC(i, width=w, flag="0"),"counts.matrix.txt",sep='')
    write.table(raw_p_data,file=raw_outfile,quote=F,sep="\t")
  }
  
  files = list.files(outdir,full.names=TRUE,include.dirs=FALSE,pattern=".matrix.txt")
  pdf(file=paste(path,prefix,file_name,"cluster_plots.pdf",sep=''))
  par(mfrow=c(2, 2))
  par(cex=0.6)
  #par(mar=c(5,3,7,2))
  for (i in 1:length(files)) {
    data = read.table(files[i], header=T, row.names=1)
    ymin = min(data); ymax = max(data);
    cluster = gregexpr("cluster_[[:digit:]]+",files[i])
    cluster = regmatches(files[i],cluster)
    plot_label = paste(prefix,cluster,', ',length(data[,1])," Features", sep='')
    par(mar = c(7, 4, 3, 0.4)+ 0.1)
    plot(as.numeric(data[1,]), type='l', ylim=c(ymin,ymax), main=plot_label, col='lightgray', 
         xaxt='n', xlab='', ylab='centered log(counts)')
    #axis(side=1, at=1:length(data[1,]), labels=colnames(data), las=2)
    axis(side=1, at=1:length(data[1,]), labels=design[,1][order(design[,1])], las=2)
    for(r in 2:length(data[,1])) {
      points(as.numeric(data[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(data)), type='o', col='blue')
  }
  dev.off()
}

getClusters1 <- function(DEdata,ldata,method=c("PAM","P","K","Km"),design,path,prefix="",w=2,val=NULL,JSD=F,hs=80){
  packages(c("clusterSim","cluster","pheatmap"))
  # Recomended options: method PAM, val=NULL
  myheatcol = colorpanel(75, 'blue','grey','red')
  data = t(scale(t(as.matrix(ldata)), scale=F)) # center rows by mean substraction
  if (JSD==T){
    gene_dist <- dist.JSD(t(DEdata))
    sample_dist = dist.JSD(DEdata)
    dis <- 'JDS'
  }else{
    gene_dist = dist(data, method='euclidean')
    sample_dist = dist(t(data), method='euclidean')
    dis <- 'Euc'
  }
  hc_genes = hclust(gene_dist, method='complete')
  sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
  sample_dist = dist(t(data), method='euclidean')
  hc_samples = hclust(sample_dist, method='complete')
  
  if (is.null(val)){
    nclusters=NULL
    
    for (k in 1:30) { 
      if (k==1) {
        nclusters[k]=NA 
      } else {
        if (method == "Km"){data.cluster_temp = kmeans(data, centers=k, iter.max=100, nstart=5)$cluster
        }else if (method == "K"){data.cluster_temp = cutree(as.hclust(hc_genes), k=k)
        }else if (method == "PAM"){data.cluster_temp = pam(as.dist(gene_dist),k, diss=TRUE)$clustering
        }
        nclusters[k]=index.G1(data,data.cluster_temp,d = gene_dist,centrotypes="medoids")
      }
    }
    pdf(paste(path,"OptClusterNum_CHindex.pdf",sep=''),width=8, height=5,onefile=FALSE)
    plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
    dev.off()
    val <- which(nclusters==max(nclusters[-1]))[1]
  }

 if (method == "Km"){
    # Use K-means clustering to define K gene sets. (use the -K parameter). 
    # This does not leverage the already hierarchically clustered genes as shown in the heatmap, 
    # and instead uses a least-sum-of-squares method to define exactly k gene clusters.
    kmeans_clustering <- kmeans(data, centers=val, iter.max=100, nstart=5)
    gene_partition_assignments = kmeans_clustering$cluster
    file_name <- paste("clusters_fixed_",dis,"_Kmeans_",val,sep="")
  }else if (method == "K"){
    # Cut the hierarchically clustered genes (as shown in the heatmap) into exactly K clusters
    gene_partition_assignments <- cutree(as.hclust(hc_genes), k=val) # cut the tree into k clusters
    file_name <- paste("clusters_fixed_",dis,"_Ktree_",val,sep="")
  }else if (method == "P"){
    #(Recommended) cut the hierarchically clustered gene tree at --Ptree percent height of the tree.
    gene_partition_assignments <- cutree(as.hclust(hc_genes), h=val/100*max(hc_genes$height))
    file_name <- paste("clusters_fixed_",dis,"_P_",val,sep="")
  }else if (method == "PAM"| is.null(method)){
  # Used by Arumugam et al 2011.Partitioning around medoids (PAM) clustering algorithm
  # to cluster the abundance profiles. PAM derives from the basic k-means algorithm, 
  # but has the advantage that it supports any arbitrary distance measure and is more robust 
  # than k-means
  gene_partition_assignments = pam(as.dist(gene_dist),val, diss=TRUE)$clustering
  file_name <- paste("clusters_fixed_",dis,"_PAM_",val,sep="")
  }
  sample_partition_assignments <- cutree(as.hclust(hc_samples), h=hs/100*max(hc_samples$height))
  file_name <- paste("sample_clusters_fixed_",dis,"_P_",hs,sep="")
  
  max_cluster_count = max(gene_partition_assignments)
  outdir <- paste(path,prefix,file_name,sep='')
  if(dir.exists(outdir)){unlink(outdir, recursive = T)}
  dir.create(outdir)
  colourCount = length(unique(gene_partition_assignments))
  #partition_colors = rainbow(colourCount, start=0.4, end=0.95)
  partition_colors <- colorRampPalette(brewer.pal(8, "Set1"))(colourCount)
  gene_colors_dframe = data.frame(clusters=gene_partition_assignments, colors=partition_colors[gene_partition_assignments])
  write.table(gene_colors_dframe, file=paste(path,prefix,file_name,".heatmap.gene_cluster_colors.dat",sep=''), quote=F, sep='  ')
  gene_colors = as.matrix(partition_colors[gene_partition_assignments])
  
  sampleColors= colorRampPalette(brewer.pal(8, "Set1"))(length(unique(sample_partition_assignments)))
  sample_colors_dframe = data.frame(clusters=sample_partition_assignments, colors=sampleColors[sample_partition_assignments])
  sample_colors = t(as.matrix(sampleColors[sample_partition_assignments]))
  
  pdf(paste(path,prefix,file_name,"unclusSamples.heatmap.pdf",sep=''))
  heatmap.3(data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=hc_samples$labels, col=myheatcol, 
            RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, margins=c(10,10),
            main=paste('Heatmap-gene clusters',sep=""),labCol=design[,1])
  dev.off()

  pdf(paste(path,prefix,file_name,".heatmap.pdf",sep=''))
  heatmap.3(data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, 
            RowSideColors=gene_colors, ColSideColors=sample_colors,scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, margins=c(10,10),
            main=paste('Heatmap-gene clusters',sep=""),labCol=design[,1])
  dev.off()
  
  gene_names = rownames(data)
  num_cols = length(data[1,])
  raw_outdir <- paste(outdir,"/rawCounts",sep='')
  if(dir.exists(raw_outdir)){unlink(raw_outdir, recursive = T)}
  dir.create(raw_outdir)
  for (i in 1:max_cluster_count) {
    partition_i = (gene_partition_assignments == i)
    partition_data = data[partition_i,,drop=F]
    #partition_data = partition_data[,hc_samples$order, drop=F]  # Uncomment if you want order samples by hcluster similarities
    partition_data = partition_data[,order(design[,1]),drop=F]
    outfile = paste(outdir,"/",prefix,"cluster_",formatC(i, width=w, flag="0"), "_log_medianCentered.matrix.txt", sep='')
    write.table(partition_data, file=outfile, quote=F, sep="\t")
    raw_p_data = DEdata[partition_i,,drop=F]
    
    #raw_p_data = raw_p_data[,raw_hc_samples$order, drop=F]      # Uncomment if you want order samples by hcluster similarities
    raw_p_data = raw_p_data[,order(design[,1]),drop=F]
    raw_outfile = paste(raw_outdir,"/",prefix,"cluster_",formatC(i, width=w, flag="0"),"counts.matrix.txt",sep='')
    write.table(raw_p_data,file=raw_outfile,quote=F,sep="\t")
  }
  
  files = list.files(outdir,full.names=TRUE,include.dirs=FALSE,pattern=".matrix.txt")
  pdf(file=paste(path,prefix,file_name,"cluster_plots.pdf",sep=''))
  par(mfrow=c(2, 2))
  par(cex=0.6)
  #par(mar=c(5,3,7,2))
  for (i in 1:length(files)) {
    data = read.table(files[i], header=T, row.names=1)
    ymin = min(data); ymax = max(data);
    cluster = gregexpr("cluster_[[:digit:]]+",files[i])
    cluster = regmatches(files[i],cluster)
    plot_label = paste(prefix,cluster,', ',length(data[,1])," Features", sep='')
    par(mar = c(7, 4, 3, 0.4)+ 0.1)
    plot(as.numeric(data[1,]), type='l', ylim=c(ymin,ymax), main=plot_label, col='lightgray', 
         xaxt='n', xlab='', ylab='centered log(counts)')
    #axis(side=1, at=1:length(data[1,]), labels=colnames(data), las=2)
    axis(side=1, at=1:length(data[1,]), labels=design[,1][order(design[,1])], las=2)
    for(r in 2:length(data[,1])) {
      points(as.numeric(data[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(data)), type='o', col='blue')
  }
  dev.off()
}

### stich from Arumugam et al 2011 - Jensen-Shannon distance:
## pulled from here, and then tweaked slightly: http://enterotype.embl.de/enterotypes.html

dist.JSD <- function(inMatrix,pseudocount=0.000001) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

net.mb <- function(df,nc=1,tl){
  packages(c("phyloseq","Matrix","igraph","devtools","RAM"))
  install_github("zdk123/SpiecEasi")
  packages(c("SpiecEasi"))
  mat <- returnAppropriateObj(df,norm=F,log=F)
  mat <- t(mat)
  a <- spiec.easi(mat,method='mb',npn=TRUE,lambda.min.ratio=1e-2,nlambda=30,icov.select.params=list(rep.num=100,ncores=nc))
  b <- graph.adjacency(a$refit, mode='undirected')
  c <- rowMeans(clr(mat, 1))+6
  d <- layout.fruchterman.reingold(b)
  #plot(b, layout=d, vertex.size=c, vertex.label=NA, main="MB")
  elist.mb <- summary(symBeta(getOptBeta(a), mode='maxabs'))
  weight <- unlist(apply(as_edgelist(b),1,function(x) {
    y <- elist.mb[elist.mb$i==x[1] & elist.mb$j==x[2],]
    if (nrow(y)!=0) return(y$x) else return(0.01)
  }))
  E(b)$weight <- abs(weight) 
  graph_attr(b,"ceb") <- cluster_edge_betweenness(b)
  graph_attr(b,"cfg") <- cluster_fast_greedy(b)
  graph_attr(b,"clp") <- cluster_label_prop(b)
  E(b)$color <- "forestgreen"
  E(b)$weight <- weight
  E(b)$color[E(b)$weight < 0] <- 'indianred1'
  graph_attr(b,"kc") <- coreness(b,mode="all")
  V(b)$Phylum <- as.vector(fData(df)$Phylum)
  levels <- taxon.level(tl,"Phylum")
  colourCount = length(levels)
  col <- colorRampPalette(brewer.pal(8, "Set1"))(9)[c(1:colourCount)]
  partition_assignments <- unlist(lapply(V(b)$Phylum, function(x) which(levels==x)))
  V(b)$pcol <- col[partition_assignments]
  graph_attr(b,"pcol") <- col
  graph_attr(b,"plev") <- levels[sort(unique(partition_assignments),decreasing=F)]
  graph_attr(b,"passort") <- assortativity_nominal(knet,partition_assignments, directed=F)
  
  V(b)$Class <- as.vector(fData(df)$Class)
  levels <- taxon.level(tl,"Class")
  colourCount = length(levels)
  palette <- c(RAM.pal(cols.needed=(10)),RAM.pal(cols.needed=(35))[c(17:21,26)],
               RAM.pal(cols.needed=(50))[c(36:50)],colorRampPalette(brewer.pal(8, "Greys")[6:1])(colourCount))
  col <- palette[c(1:colourCount)]
  partition_assignments <- unlist(lapply(V(b)$Class, function(x) which(levels==x)))
  V(b)$ccol <- col[partition_assignments]
  graph_attr(b,"clev") <- levels[sort(unique(partition_assignments),decreasing=F)]
  graph_attr(b,"ccol") <- col
  graph_attr(b,"cassort") <- assortativity_nominal(knet,partition_assignments, directed=F)
  
  V(b)$Family <- as.vector(fData(df)$Family)
  levels <- taxon.level(tl,"Family")
  colourCount = length(levels)
  palette <- c(RAM.pal(cols.needed=(10)),RAM.pal(cols.needed=(35))[c(17:21,26)],
               RAM.pal(cols.needed=(50))[c(36:50)],colorRampPalette(brewer.pal(8, "Greys")[6:1])(colourCount))
  col <- palette[c(1:colourCount)]
  partition_assignments <- unlist(lapply(V(b)$Family, function(x) which(levels==x)))
  V(b)$fcol <- col[partition_assignments]
  graph_attr(b,"flev") <- levels[sort(unique(partition_assignments),decreasing=F)]
  graph_attr(b,"fcol") <- col
  graph_attr(b,"fassort") <- assortativity_nominal(knet,partition_assignments, directed=F)
  
  V(b)$Genus <- as.vector(fData(df)$Genus)
  levels <- taxon.level(tl,"Genus")
  colourCount = length(levels)
  palette <- c(RAM.pal(cols.needed=(10)),RAM.pal(cols.needed=(35))[c(17:21,26)],
               RAM.pal(cols.needed=(50))[c(36:50)],colorRampPalette(brewer.pal(8, "Greys")[6:1])(colourCount))
  col <- palette[c(1:colourCount)]
  partition_assignments <- unlist(lapply(V(b)$Genus, function(x) which(levels==x)))
  V(b)$gcol <- col[partition_assignments]
  graph_attr(b,"glev") <- levels[sort(unique(partition_assignments),decreasing=F)]
  graph_attr(b,"gcol") <- col
  graph_attr(b,"gassort") <- assortativity_nominal(knet,partition_assignments, directed=F)
  
  return(list(mb=a,graph=b,layout=d,vsize=c))
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

twographs <- function(){
  library(ggplot2)
  library(gridExtra)
  p <- "/Users/guillermotorres/Documents/OtrosUsuarios/Lena/" ## path to the folder where the input files are.
  gene.p <- read.table(paste(p,"genes_with_percentages.txt",sep=''),sep="\t")
  gene.p$P <- as.numeric(gsub('%','',gene.p$V2,fixed=T))
  gene.p <- gene.p[order(-gene.p$P),]
  gene.p$V1 <- factor(as.character(gene.p$V1),levels=as.character(gene.p$V1))
  gene.names <- unlist(apply(gene.p,1,function(x) paste(x[1],x[2],collapse="  ")))
  names(gene.names) <- levels(gene.p$V1)
  
  ind.load <- read.table(paste(p,"mutational_load_per_sample.100spBCCs.txt",sep=''),sep="\t")
  ind.load <- ind.load[order(-ind.load$V2),]
  ind.load$V1 <- factor(as.character(ind.load$V1),levels=as.character(ind.load$V1))
  
  ind.gen <- read.table(paste(p,"merged_data.genes_fig1_+_candidates.txt",sep=''),sep="\t")
  ind.gen$V2 <- factor(ind.gen$V2,levels=levels(gene.p$V1))
  ind.gen$V1 <- factor(ind.gen$V1,levels=levels(ind.load$V1))
  
  g1 <- ggplot(ind.load,aes(x=V1,y=V2))+geom_bar(stat='identity')+
    scale_y_discrete("Mutational load",limits=c(seq(0,max(ind.load$V2),length.out=5)))+
    theme(axis.line.x = element_blank(),
          axis.line.y = element_line(color="black",size=0.5),
          axis.text.x=element_blank(),axis.title.x=element_blank(),#axis.text.x= element_text(angle = 90, hjust = 1),
          panel.background=element_blank(),legend.position="none",axis.title.y=element_text(margin=margin(0,20,0,0)))
  
  g2 <- ggplot(ind.gen,aes(x=V1,y=V2,fill=V3))+geom_point(shape=22)+
    scale_y_discrete("Genetic alteration",labels=gene.names)+
    facet_grid(V2 ~ ., scales = "free", space = "free") +
    theme(axis.line.x = element_blank(),
          axis.text.x=element_blank(),axis.title.x=element_blank(),#axis.text.x= element_text(angle = 90, hjust = 1),#
          panel.background=element_blank(),
          #panel.grid.major.y = element_line(size = 3,color='white',linetype=1),
          panel.grid.major.x = element_line(size=2,color='grey88',linetype=1),
          legend.position="none",strip.background = element_blank(),
          strip.text = element_blank(),axis.title.y=element_text(margin=margin(0,20,0,0)))
  ## legend purposes ##
  gl <- ggplot(ind.gen,aes(x=V1,y=V2,fill=V3))+geom_point(shape=22)+
    scale_y_discrete("Genetic alteration",labels=gene.names)+
    facet_grid(V2 ~ ., scales = "free", space = "free") +
    theme(axis.line.x = element_blank(),
          axis.text.x=element_blank(),axis.title.x=element_blank(),#axis.text.x= element_text(angle = 90, hjust = 1),#
          panel.background=element_blank(),
          #panel.grid.major.y = element_line(size = 3,color='white',linetype=1),
          panel.grid.major.x = element_line(size=2,color='grey88',linetype=1),
          strip.background = element_blank(),
          strip.text = element_blank(),legend.title=element_blank())
  ##
  gp1 <- ggplot_gtable(ggplot_build(g1))
  gp2 <- ggplot_gtable(ggplot_build(g2))
  gpl <- ggplot_gtable(ggplot_build(gl))
  gp.leg <- gpl$grobs[[56]]               ## extracting legend information - number depends on the amount of elements in the graph
  maxWidth = grid::unit.pmax(gp1$widths[2:3],gp2$widths[2:3])
  gp1$widths[2:3] <- as.list(maxWidth)
  gp2$widths[2:3] <- as.list(maxWidth)
  graph <- grid.arrange(arrangeGrob(gp1,gp2,heights=c(2/6,4/6),ncol=1),gp.leg,widths=c(9/10,1/10))
  #graph ## uncoment this line to visualize the plot in r
  ggsave(paste(p,'myplot.pdf',sep=''),graph,width=12, height=8,device="pdf")
  
  
  
  gp.leg <- gp2$grobs[[8]]
  
  hist <- data.frame(id=c(1:10),mtmo=sample(10:100,10,replace=T))
  ids <- as.factor(hist$id)
  hist$id <- ids
  genmut <- data.frame(id=c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,5,5,5,6,6,6,6,7,7,8,9,10),
              gen=c('a','b','c','a','b','c','d','e','a','c','d','a','c','d','e','a','b','c',
                    'b','c','d','e','a','b','c','d','e'),
              mut=c(rep('T',10),'Non','Non','Non',rep('T',10),'Non','Non','T','T'))
  genmut$id <- factor(genmut$id,levels=levels(ids))
  hist$table <- rep('t1',nrow(hist))
  genmut$table <- rep('t2',nrow(genmut))
  
  g1 <- ggplot(hist,aes(x=as.factor(id),y=mtmo))+geom_bar(stat='identity')+
    scale_y_discrete(limits=c(seq(0,100,length.out=5),120))+
    theme(axis.line.x = element_blank(),
          axis.line.y = element_line(color="black",size=0.5),
          axis.text.x= element_blank(),axis.title.x=element_blank(),
          panel.background=element_blank(),legend.position="none")
  g2 <- ggplot(genmut,aes(x=id,y=gen,fill=mut))+geom_point(shape=22, size=10)+
    theme(axis.line.x = element_blank(),
          axis.line.y = element_line(color="black",size=0.5),
          axis.text.x= element_blank(),axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.grid.major.x = element_line(size = 1,color='grey'),
          legend.position="none")
  gp1 <- ggplot_gtable(ggplot_build(g1))
  gp2 <- ggplot_gtable(ggplot_build(g2))
  gp.leg <- gp2$grobs[[8]]
  maxWidth = grid::unit.pmax(gp1$widths[2:3],gp2$widths[2:3])
  gp1$widths[2:3] <- as.list(maxWidth)
  gp2$widths[2:3] <- as.list(maxWidth)
  grid.arrange(arrangeGrob(gp1,gp2,heights=c(2/6,4/6),ncol=1),gp.leg,widths=c(9/10,1/10))
  
  ## http://stackoverflow.com/questions/15999304/plotting-continuous-and-discrete-series-in-ggplot-with-facet
  
  
  ggplot()+
    geom_bar(data=hist,aes(x=id,y=mtmo),stat='identity')+
    #scale_y_discrete(data=hist,limits=c(seq(0,100,length.out=5),120))+
    facet_grid(table~.,scale='free_y')+
    geom_point(data=genmut,aes(x=id,y=gen,fill=mut),shape=22, size=10)+
    theme(axis.line.x = element_blank(),
          axis.line.y = element_line(color="black",size=0.5),
          axis.text.x= element_blank(),axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.grid.major.x = element_line(size = 1,color='grey'))
  
}

#variable <- 'Salinity'#opt$variable
#conf <- c('Arena','Fe')#opt$conf
#pval <- 0.5 # opt$pval
#out <- paste(ra,'humann2/predicted_results/',sep='') #opt$out
#level <- 'modulesXMP' #opt$level
#level <- levelm
#df <- dfk
fitzigdiff <- function(df,variable,conf,pval,out,level){

  ## Differential abbundance test - including counfounders ###
  ## FitZig Method
  p.f <- cumNormStat(df,pFlag=TRUE,main="Data") # Calculates the percentile for which to sum counts up to and scale by.
  df <- cumNorm(df,p=p.f)
  try(cfs <- unlist(sapply(conf,function(x) return(paste('pData(df)$',x,sep='')))))# Confounders
  
  nf <- normFactors(df)  # Calculates each column's quantile and calculates the sum up to and including p quantile
  normFactor <- normFactors(df)#
  if (length(conf)==1 & conf==''){
    mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df)[["',variable,'"]])+',sep=''),'normFactor')))
  }else{
    mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df)[["',variable,'"]])+',sep=''),
                                         paste(cfs,collapse='+'),'+normFactor')))
  }
  mod <- model.matrix(as.formula(paste("~",paste('0+factor(pData(df)[["',variable,'"]])+',sep=''),
                                       paste(cfs,collapse='+'))))
  settings <- zigControl(maxit=10,verbose=T)
  fit <- NULL
  try(fit <- fitZig(obj=df,mod=mod,control=settings,useMixedModel=T,useCSSoffset=T))
  if(length(fit)==0){
    fit <- fitZig(obj=df,mod=mod,control=settings,useMixedModel=F)
  }
  #summary(calculateEffectiveSamples(fit))
  zigFit <- fit$fit
  finalMod <- fit$fit$design
  x <- unlist(lapply(colnames(finalMod),function(x) return(gsub(paste('factor(pData(df)[["',variable,'"]])',sep=''),variable,x,fixed=T))))
  colnames(finalMod) <- x
  x <- unlist(lapply(colnames(finalMod),function(x) return(gsub('pData(df)$','',x,fixed=T))))
  colnames(finalMod) <- x
  colnames(zigFit$coefficients) <- colnames(finalMod)
  colnames(zigFit$stdev.unscaled) <- colnames(finalMod)
  colnames(zigFit$cov.coefficients) <- colnames(finalMod)
  colnames(zigFit$design) <- colnames(finalMod)
  k <- levels(factor(pData(df)[[variable]]))  # states of the variable in study
  k <- unlist(sapply(k,function(x) return(paste(variable,x,sep=''))))
  mk <- combn(k,2)
  contrast_list <- c()
  for (i in 1:dim(mk)[2]){contrast_list <- c(contrast_list,paste(mk[1,i],mk[2,i],sep='-'))}
  contrast.matrix <- makeContrasts(contrasts=contrast_list,levels=finalMod)
  fit2 <- contrasts.fit(zigFit,contrast.matrix)
  fit2=eBayes(fit2)
  results <- decideTests(fit2,method="separate",adjust.method="fdr",p.value=pval)#,lfc=0.01) ## 1/0/-1 result DEtable 
  DElist <- c()
  for (i in 1:length(colnames(fit2$coef))){
    lm <- rownames(subset(results,results[,i]!=0))
    if (length(lm)>0){
      message (colnames(fit2$coef)[i]," count talbes generation...")
      texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=pval,adjust="BH")
      if (nrow(texp)>0){
        texp$DE <- results[rownames(results)%in%rownames(texp),i]
        #texp <- cbind(texp,fData(df)[rownames(texp),-c(1,2)])
        try(texp <- cbind(texp,fData(df)[rownames(texp),2]))
        ipath <- data.frame(pathway=rownames(texp),colors=ifelse(texp$DE==1,'#FF0000','#0000FF'),width='W2')
        write.table(ipath,file=paste(out,level,'vs',colnames(fit2$coef)[i],'_conf_',paste(conf,collapse='_'),'_fitZig_ipath.txt',sep=''),
                    quote=F,sep="\t",na="NA",row.names=F)
        write.table(texp,file=paste(out,level,'vs',colnames(fit2$coef)[i],'_conf_',paste(conf,collapse='_'),'_fitZig.txt',sep=''),
                    quote=F,sep="\t",na="NA",row.names=T)
        DElist <- c(DElist,rownames(texp))
      }
    }else{
      message(paste(colnames(fit2$coef)[i],' has no significant elements at ',pval,' threshold'))
      #texp <- topTable(fit2,coef=i,number=NROW(fit2$coefficients),p.value=opt$pval,adjust="BH")
      #write.table(texp,file=paste(opt$out,opt$level,'vs',colnames(fit2$coef)[i],'_conf_',paste(opt$conf,collapse='_'),'_fitZig.txt',sep=''),
      #            quote=F,sep="\t",na="NA",row.names=T)
    }
  }
  DElist.1 <- unique(DElist)
  f <- which(rownames(MRcounts(df))%in%DElist.1)
  dfc <- df[f,1:length(sampleNames(df))]
  dfexport_counts <- setNames(data.frame(round(MRcounts(dfc,norm=T),2)),colnames(MRcounts(dfc)))
  write.table(dfexport_counts,file=paste(out,level,'DifAbund_conf_',paste(opt$conf,collapse='_'),'_Counts.txt',sep=''),
                          quote=F,sep="\t",na="NA",row.names=T)
  #try(dfexport_alltax <- fData(df)[rownames(dfexport_counts),2,drop=F])
  #if(level=='OTU'){
  #  dfexport_alltax <- lapply(rownames(dfexport_counts),function(x){
  #    fData(df.r)[fData(df.r)['OTU']==x,][1,3:8]
  #  })
  #}else{
  #  dfexport_alltax <- lapply(rownames(dfexport_counts),function(x){
  #    fData(df.r)[fData(df.r)[opt$level]==x,][1,3:which(colnames(fData(df.r))==opt$level)]
  #  })
  #}
#  
  #dfexport_alltax.df <- do.call(rbind.data.frame,dfexport_alltax)
  #dfexport <- cbind(dfexport_counts,dfexport_alltax)
  #write.table(dfexport,file=paste(out,level,'DifAbund_conf_',paste(opt$conf,collapse='_'),'_CountsTaxonomy.txt',sep=''),
  #            quote=F,sep="\t",na="NA",row.names=T)
  
  #### Cluster analysis
  ##source(toolbox)
  #if(nrow(MRcounts(dfc))>0){
  #  if(!is.na(refvar)) {design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,refvar])),c('ID',opt$variable))
  #  }else design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,opt$variable])),c('ID',opt$variable))
  #  design <- design[with(design,order(design[,2])),]
  #  rownames(design) <- design$ID
  #  getClusters( MRcounts(dfc,norm=T), MRcounts(dfc,norm=T,log=T),
  #               design=design,method=opt$clmethod,path=opt$out,
  #               prefix=paste(opt$level,'Sign_conf_',paste(opt$conf,collapse='_'),'_',sep=""),val=opt$clval,
  #               cellwidth=20,cellheight=8.3,text_size=8.7,height=13,width=7)
  #}
  ### ploting all genders..
  ##if(!is.na(refvar)) {design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,refvar])),c('ID',opt$variable))
  ##}else design <- setNames(data.frame(rownames(pData(dfc)),as.factor(pData(dfc)[,opt$variable])),c('ID',opt$variable))
  ##design <- design[with(design,order(design[,2])),]
  ##rownames(design) <- design$ID
  ##getClusters(MRcounts(df,norm=T), MRcounts(df,norm=T,log=T),
  ##             design=design,method=opt$clmethod,path=opt$out,prefix=paste(opt$level,'All_',sep=""),val=50)
#  
  #} 

}




###
###
# Geting sourcing paths
#initial.options <- commandArgs(trailingOnly = FALSE)
#file.arg.name <- "--file="
#script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
#script.basename <- dirname(script.name)
#other.name <- paste(sep="/", script.basename, "toolbox.R")
#print(paste("Sourcing",other.name,"from",script.name))
###

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# Borrowing the following from gplots  (gplots isn't compatible with R 3.0 (yet), and so bypassing it for now).

colorpanel = function (n, low, mid, high) 
{
  if (missing(mid) || missing(high)) {
    low <- col2rgb(low)
    if (missing(high)) 
      high <- col2rgb(mid)
    else high <- col2rgb(high)
    red <- seq(low[1, 1], high[1, 1], length = n)/255
    green <- seq(low[3, 1], high[3, 1], length = n)/255
    blue <- seq(low[2, 1], high[2, 1], length = n)/255
  }
  else {
    isodd <- odd(n)
    if (isodd) {
      n <- n + 1
    }
    low <- col2rgb(low)
    mid <- col2rgb(mid)
    high <- col2rgb(high)
    lower <- floor(n/2)
    upper <- n - lower
    red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
                                                                1], high[1, 1], length = upper))/255
    green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
               seq(mid[3, 1], high[3, 1], length = upper))/255
    blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
              seq(mid[2, 1], high[2, 1], length = upper))/255
    if (isodd) {
      red <- red[-(lower + 1)]
      green <- green[-(lower + 1)]
      blue <- blue[-(lower + 1)]
    }
  }
  rgb(red, blue, green)
}


greenred = function (n)  {
  colorpanel(n, "green", "black", "red")
}

odd = function (x) {
  x%%2 == 1
}

even = function (x) {
  x%%2 == 0
}


heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.1,
                      #cexRow = 0.2 + 1/log10(max(nr,2)),
                      #cexCol = 0.2 + 1/log10(max(nc,2)),
                      cexRow = 0.2,
                      cexCol = 0.2,                  
                      
                      scaleRangeMin,
                      scaleRangeMax,
                      
                      
                      cex.main = 1,
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  retval <- list()
  
  
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  
  dendrogram <- match.arg(dendrogram)
  
  trace <- match.arg(trace)
  
  density.info <- match.arg(density.info)
  
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  
  nr <- di[1]
  nc <- di[2]
  
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  #print(paste("nr:", nr, "nc:", nc, "cexCol:", cexCol, "cexRow:", cexRow))
  #stop("debug")
  
  
  
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  
  x <- x[rowInd, colInd]  # rearrange matrix according to dendrograms
  x.unscaled <- x
  
  cellnote <- cellnote[rowInd, colInd]  # also rearrange the cellnotes
  
  # get labels 
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  
  
  ## do scaling of matrix according to Z-scores
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  
  # number of breaks
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  
  # set breakpoints
  if (length(breaks) == 1) {
    if (missing(scaleRangeMin))
      scaleRangeMin = min(x, na.rm=na.rm)
    
    if (missing(scaleRangeMax))
      scaleRangeMax = max(x, na.rm=na.rm)
    
    
    if (!symbreaks) {
      breaks <- seq(scaleRangeMin, scaleRangeMax, length=breaks);
    } else {
      #extreme <- max(abs(x), na.rm = TRUE)
      extreme = max(abs(c(scaleRangeMin,scaleRangeMax)), na.rm=na.rm)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  
  if (class(col) == "function")
    col <- col(ncol)
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  
  # adjust for out-of-range given break settings
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  
  # layout height
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  
  # layout width
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  
  # define the layout
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || ncol(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of ncol(x) ", nc, " columns")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
      side_height = min(side.height.fraction*nrow(ColSideColors), 1);
      lhei=c(lhei[1], side_height, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || nrow(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of nrow(x) ", nr, " rows.  It currently has ", nrow(RowSideColors), " rows.")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
      side_width = min(side.height.fraction*ncol(RowSideColors), 1);
      lwid <- c(lwid[1], side_width, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  ###########################################
  ## Draw the colorbars for the annotations:
  ###########################################	
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[rowInd, , drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      # print(rsc)
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      #print("RSC: ", rsc)
      #print(rsc.colors)    
      image(1:nrow(rsc), 1:ncol(rsc), rsc, col = as.vector(rsc.colors), axes = FALSE, xlab="", ylab="")
      
      # add labels
      if (length(colnames(RowSideColors)) > 0) {  
        #axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
        #axis(1, 0:(nrow(rsc)-1), colnames(RowSideColors), las = 2, tick = T) # ncol because transposed
        axis(1, 1:ncol(RowSideColors), labels=colnames(RowSideColors), las=2, cex.axis=0.5, tick=F, xlab="", ylab="")
        
      }
    }
  }
  
  
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[, colInd, drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      #print(csc)
      image(1:nrow(t(csc)), 1:ncol(t(csc)), t(csc), col = as.vector(csc.colors), axes = FALSE, xlab="", ylab="")
      
      # add labels
      if (length(rownames(ColSideColors)) > 0) {
        #axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        axis(2, 1:(nrow(ColSideColors)), labels=rownames(ColSideColors), las = 2, tick = FALSE, cex.axis=0.5)
      }
    }
  }
  
  
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  
  # draw the central heatmap
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  
  # store the matrix drawn
  retval$carpet <- x
  
  # store the dendrograms
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  
  # store the breaks
  retval$breaks <- breaks
  
  # store the colormap used
  retval$col <- col
  
  # specially color in the na values	
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color, add = TRUE)
  }
  
  # X-axis column labels
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
  
  # X-axis title
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  
  # Y-axis row labeling
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  
  # Y-axis title
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  
  
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  
  # column trace
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  
  # row trace
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  
  # add cell labels
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)
  
  ###########################
  ## Plot the row dendrogram
  ###########################
  
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  
  #############################
  ## Plot the column dendrogram
  #############################
  
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  
  if (!is.null(main))
    title(main, cex.main=cex.main) #cex.main = 1.5 * op[["cex.main"]])
  
  
  ############################
  ## Add the Color Chart
  ############################
  
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(c(x,breaks), na.rm = TRUE)
      max.raw <- max(c(x,breaks), na.rm = TRUE)
    }
    
    message('for plotting:: min.raw: ', min.raw, ' max.raw: ', max.raw);
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], high = retval$breaks[-1], color = retval$col)
  
  invisible(retval)
}


