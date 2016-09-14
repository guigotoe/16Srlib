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
    setRepositories(ind=1:10)
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has])
  }
  lapply(requirements, require, character.only = TRUE)
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
  df[[col]] <- gsub(";unclassified",'',df[[col]])
  df[[col]] <- gsub(";",'',df[[col]])
  return(df)
}

mothur.taxonomy <- function(taxonomy){
  packages(c("reshape2"))
  taxnames <- c('Kingdom','Phylum','Class','Order','Family','Genus',"Specie","others")
  x <- colsplit(taxonomy$Taxonomy,pattern="\\([[:digit:]]*\\);",names = taxnames)
  for (i in taxnames){
    x <- changename(x,i)
  }
  x[x==""|is.na(x)]  <- 'unclassified'
  df <- cbind(taxonomy[,1:2],x)
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
  pdf(paste(o,"Heatmap_",gsub(" ",'_',title),".pdf",sep=""),pointsize=12, width=10,height=7)
  pheatmap(sampleDistMatrix,main=title,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
}

taxprop <- function(dfp,v,tl,o,limit=0.01,u=F){
  ## input: 1) MRexperiment (proportions,features,metadata) 2) reference variable: x axis 3) taxlevel 4) output path
  dfx <- MRcounts(dfp)
  tx <- fData(dfp)
  mx <- pData(dfp)
  rownames(dfx) <- tx[[tl]]
  dfxt<- as.data.frame(t(dfx))
  ## sum columns with same name -- OTUs with same taxonomical level:
  dfxta <- as.data.frame(sapply(unique(names(dfxt)[duplicated(names(dfxt))]), function(x) Reduce("+", dfxt[ , grep(x, names(dfxt))])))
  gm <- cbind(dfxta,mx[v],id=mx[,1])
  gm_m <- melt(gm,id.vars=c(v,"id"))
  taxa <- unlist(apply(gm_m,1, function(x) if(as.numeric(x["value"])>limit){(x["variable"])}else{"Others"}))
  gm_m$"taxa" <- as.factor(taxa)
  gmx2 <- aggregate(as.formula(paste("value~",'id+taxa+',paste(v,collapse="+"))),data=gm_m,FUN=sum) # Because some taxa now are xOthers so we need to summ their values
  
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
  
  title <- paste(tl," abundance distribution",sep='')
  write.table(gmx2,file=paste(o,gsub(" ",'_',title),'_',v,'.txt',sep=''),sep="\t",quote=F,row.names=F)
  ut <- ''
  if(u==T){gmx2 <- subset(gmx2,gmx2$taxa!='unclassified');ut <- '_Xuncl'}
  
  distplot <- ggplot(gmx2,aes(x=as.factor(gmx2[[v]]),y=value,fill=taxa))+geom_boxplot()+#geom_jitter(position=position_jitter(width=.2), size=0.5)+
    scale_fill_manual(name=tl,values=palette) +
    scale_x_discrete(v)+ylab("Proportion")+
    ggtitle(paste(title,sep=""))+
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
          panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          legend.position="bottom",legend.box="horizontal",
          legend.text = element_text(size=8),
          legend.title = element_text(size=10, face="bold"),
          plot.title = element_text(lineheight=.12, face="bold"))
  ggsave(paste(gsub(" ",'_',title),ut,'_',v,'.pdf',sep=''),plot=distplot,path=o,width=8,height=5,device="pdf")
  
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

