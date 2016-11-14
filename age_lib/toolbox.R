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
#DEdata <- counts
#ldata <- counts.l
#method <- opt$clmethod
#design <- pData(df.f)[c("Age",opt$variable)]
#path <- opt$out
#prefix <- paste(opt$level,'_',sep="")
#val=opt$clval


getClusters <- function(DEdata,ldata,method=c("PAM","P","K","Km"),design,path,prefix="",w=2,val=NULL,JSD=F,hs=80){
  packages(c("clusterSim","cluster"))
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

net.mb <- function(df){
  packages(c("phyloseq","Matrix","igraph","SpiecEasi"))
  #install_github("zdk123/SpiecEasi")
  mat <- returnAppropriateObj(df,norm=T,log=F)
  mat <- t(mat)
  a <- spiec.easi(mat,method='mb',npn=TRUE,lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
  b <- graph.adjacency(a$refit, mode='undirected')
  c <- rowMeans(clr(mat, 1))+4
  d <- layout.fruchterman.reingold(b)
  #plot(b, layout=d, vertex.size=c, vertex.label=NA, main="MB")
  elist.mb <- summary(symBeta(getOptBeta(a), mode='maxabs'))
  weight <- unlist(apply(as_edgelist(b),1,function(x) {
    y <- elist.mb[elist.mb$i==x[1] & elist.mb$j==x[2],]
    if (nrow(y)!=0) return(y$x) else return(0.01)
  }))
  E(b)$weight <- weight
  E(b)$color <- "blue"
  E(b)$color[E(b)$weight < 0] <- 'red'
  return(list(mb=a,graph=b,layout=d,vsize=c))
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


