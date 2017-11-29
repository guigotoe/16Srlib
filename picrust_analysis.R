####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: November 2017
# Created: November 2017
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript picrust_analysis.R 
#
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
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/age_lib/toolbox.R"
source(toolbox)

packages(c("metagenomeSeq","optparse"))

## Options ##
#pa = '/home/torres/ikmb_storage/Metagenome/16s/2017/4_2017/'                   # AGING VERSION
#ra <- '/home/torres/Documents/Projects/Metagenome/2017_results/9_2017/'        # AGING VERSION
pa <- '/home/torres/ikmb_storage/Mangrove/16Sfa/2016/mothur_gg/'                   # MANGROVE VERSION
ra <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/'          # MANGROVE VERSION

option_list <- list(
  make_option(c("-i","--input"),action="store",type="character",default=paste(pa,'16s_OTU_table.biom',sep=''),#'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared'
              help="Path to input biom input file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(pa,'metadata.txt',sep=''),
              help="Path to input metadata file"),
  make_option(c("-c","--contribution"),action="store",type="character",default=paste(ra,'ko_contributions.txt',sep=''),
              help="Path to contribution file"),
  make_option(c("-d","--diff"),action="store",type="character",default=paste(ra,'GenusDifAbund_conf__CountsTaxonomy.txt',sep=''),
              help="Path to differentially abundant Genus file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(ra,'',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-t","--shared"),type="double",default=0.2,
              help="OTU presence; percentage of samples sharing each OTU. 0-1; default: %default"),
  make_option(c("-d","--depth"),type="integer",default=1,
              help="Minimum depth count; default: %default counts per otu"),
  make_option(c("-L","--libraries"),type="integer",default=3,
              help="Minimum depth count per libraries [depth,j]: cutoff at 'depth' counts in at least j libraries."),
  make_option(c("-r","--replicates"),action="store_true",default=FALSE,
              help="Technical replications")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$input)){stop(sprintf("There is not counts biom-file specified"))
}else if(is.na(opt$metadata)){stop(sprintf("There is not metadata file specified"))}
if(length(grep("/$",opt$out))==0) opt$out <- paste(opt$out,"/",sep="")

picrust_path <- '/home/torres/Bin/picrust-1.1.2/scripts/' # local
picrust_outpath <- paste(opt$out,'/picrust/analisys/',sep='')

if(dir.exists(picrust_outpath)){message('picrustOut folder already exist, files will be overwritten')
}else dir.create(picrust_outpath,showWarnings=F)

###### end ######
packages(c('biomformat','ggplot2','dplyr'))
ggdata <- mothur.biom(opt$input,opt$metadata)
pData(ggdata) <- pData(ggdata)[,colSums(is.na(pData(ggdata)))==0] ## remove coluns with NAs
pData(ggdata)$ID_ref <- factor(as.character(pData(ggdata)$ID_ref),levels=c('low','med','high'))
p.f <- cumNormStatFast(ggdata) # Calculates the percentile for which to sum counts up to and scale by.
ggdata <- cumNorm(ggdata,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(ggdata)
head(fData(ggdata))

## humann ####
predicted.r <- read.table(paste(ra,'picrust/ko_predicted_metagenomes.txt',sep=''),header=T,sep="\t",quote="",check.names=F)
write.table(predicted.r[,-ncol(predicted.r)],paste(ra,'humann2/ko_predicted.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)


#salloc --job-name=microb --nodes=1 --cpus-per-task=4 --partition=long --mem=16000 
# mkdir split_files
# humann2_split_table -i ko_biom_predicted_metagenomes.biom -o split_files/ --taxonomy_index -1
#humann2_split_table -i ko_contributions.txt -o split_files/ --taxonomy_index -1

# mkdir outL

#for tsv in ./split_files/*.tsv; do humann2 -i $tsv -o outL/ --pathways-database /home/sukmb347/bin/humann/data/keggc; done  #--minpath off

#humann2_join_tables -i ./out_contrib/ -o ./humann2_pathabundance.txt --file_name pathabundance
#humann2_join_tables -i ./out_contrib/ -o ./humann2_pathcoverage.txt --file_name pathcoverage

# for biom in split_files/*.biom; do humann2 -i $biom -o outL/ --pathways-database /home/sukmb347/bin/humann/data/keggc; done

# humann2_join_tables -i out/ -o ./humann2_pathabundance.txt --file_name pathabundance
# humann2_join_tables -i out/ -o ./humann2_pathcoverage.txt --file_name pathcoverage
# sed 's/_Abundance//g' humann2_pathabundance.txt > humann2_pathabundance_fixed.txt
# sed 's/_Coverage//g' humann2_pathabundance.txt > humann2_pathabundance_fixed.txt
# normalizing
# humann2_renorm_table -i humann2_pathabundance.txt -o humann2_pathabundance_relab.txt -n relab
# sed 's/_Abundance//g' humann2_pathabundance_relab.txt > humann2_pathabundance_relab_fixed.txt
# humann2_renorm_table -i humann2_pathabundance.txt -o humann2_pathabundance_cpm.txt -n cpm
# sed 's/_Abundance//g' humann2_pathabundance_cpm.txt > humann2_pathabundance_cpm_fixed.txt
### end humann##

## humann2 ##

sample_ref <- pData(ggdata)
sample_ref[['ID_ref']] <- factor(as.character(pData(ggdata)[['ID_ref']]),levels=c('low','med','high'))

## presence/absense -- get modules with presence in the samples 
mopa.r <- read.table(paste(ra,'humann2/predicted_results/Coverage_predicted_modulec_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
head(mopa.r)
mopa <- mopa.r[!mopa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),2:ncol(mopa.r)]
rownames(mopa) <- as.character(mopa.r[!mopa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),1])
## filtrating only those ko by ID_ref mean is > 0.5
th <- 0.5
sigmopa <- c()
for (i in levels(sample_ref[['ID_ref']])){
  sample <- apply(as.matrix(mopa[,as.character(sample_ref[sample_ref[['ID_ref']]==i,1])]),2,as.numeric)
  sigmopa <- c(sigmopa,rownames(mopa)[which(rowMeans2(sample)>th)])
}  
sigmopa <- unique(sigmopa)
write.table(sigmopa,paste(ra,'humann2/predicted_results/sigMO_pa.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)

## modules abundance --- just evaluate DIFF abundant with positive presence (prev. step). 
moa.r <- read.table(paste(ra,'humann2/predicted_results/Abundance_predicted_modulec_relab_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
tail(moa.r)
moa <- moa.r[!moa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),2:ncol(moa.r)]
rownames(moa) <- as.character(moa.r[!moa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),1])
#moa <- moa.r[,2:ncol(mopa.r)]
#rownames(moa) <- as.character(moa.r[,1])
# 
mo_info <- read.table(paste('/home/torres/Bin/humann/data/map_kegg.txt',sep=''),header=F,sep="\t",quote='')
mo_data <- mo_info[mo_info$V1%in%rownames(moa),]
#mo_data <- rbind(mo_data,data.frame(row.names=c('UNINTEGRATED','UNMAPPED'),V1=c('UNINTEGRATED','UNMAPPED'),V2=rep(NA,2)))
rownames(mo_data) <- mo_data[,1]
df <- newMRexperiment(decostand(moa,method='hellinger'),AnnotatedDataFrame(sample_ref),AnnotatedDataFrame(mo_data))
#df <- newMRexperiment(moa,AnnotatedDataFrame(sample_ref),AnnotatedDataFrame(mo_data))
variable <- 'Salinity'#opt$variable
conf <- c('Arena')#opt$conf
pval <- 0.05 # opt$pval
out <- paste(ra,'humann2/predicted_results/',sep='') #opt$out
level <- 'modulesXMP' #opt$level
fitzigdiff(df,variable,conf,pval,out,level)



koa.r <- read.table(paste(ra,'humann2/predicted_results/Abundance_predicted_keggc_cpm_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
head(koa.r)
koa <- koa.r[!koa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),2:ncol(koa.r)]
rownames(koa) <- gsub('ko','',as.character(koa.r[!koa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),1]),fixed=T)
#moa <- moa.r[,2:ncol(mopa.r)]
#rownames(moa) <- as.character(moa.r[,1])
# 
ko_info <- read.csv(paste('/home/torres/Bin/humann/data/ecc',sep=''),header=F,sep="\t",quote='',as.is=T)
#ko_data <- ko_info[ko_info$V1%in%rownames(koa),]
#mo_data <- rbind(mo_data,data.frame(row.names=c('UNINTEGRATED','UNMAPPED'),V1=c('UNINTEGRATED','UNMAPPED'),V2=rep(NA,2)))
#rownames(ko_data) <- ko_data[,1]

dfk <- newMRexperiment(decostand(koa,method='hellinger'),AnnotatedDataFrame(sample_ref))
#dfk <- newMRexperiment(koa,AnnotatedDataFrame(sample_ref))
write.table(MRcounts(dfk),paste(ra,'humann2/mapXMP.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)
levelm <- 'mapXMP'
conf <- c('Arena')
fitzigdiff(dfk,variable,conf,pval,out,levelm)

tail(predicted.r)
predicted <- apply(predicted.r[,-c(1,ncol(predicted.r))],2,as.numeric)
head(predicted)
cnames <- as.character(predicted.r[-which(rowSums2(predicted)==0),1])
predicted <- predicted[-which(rowSums2(predicted)==0),]
rownames(predicted) <- cnames
dfp <- newMRexperiment(decostand(predicted,method='hellinger'),AnnotatedDataFrame(sample_ref))
#dfp <- newMRexperiment(predicted,AnnotatedDataFrame(sample_ref))
levelk <- 'KO'
fitzigdiff(dfp,variable,conf,pval,out,levelk)

source(toolbox)

kopa.r <- read.table(paste(ra,'humann2/predicted_results/Abundance_predicted_keggc_relab_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
head(kopa.r)
#colnames(kocov.r) <- gsub('_Coverage','',colnames(coverage.r),fixed=T)
#tokeep.coverage <- unlist(lapply(as.character(kocov.r[,1]),function(x) length(grep('|',x,fixed=T))==0))
#kocov.r.pa <- kocov.r[tokeep.coverage,]
tpdata <- as.data.frame(t(pData(ggdata)),stringsAsFactors=FALSE)
tpdata$Pathway <- as.character(rownames(tpdata))
tpdata <- tpdata[,c(ncol(tpdata),1:(ncol(tpdata)-1))]
kopa <- kopa.r[!kopa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),2:ncol(kopa.r)]
rownames(kopa) <- as.character(kopa.r[!kopa.r$`# Pathway`%in%c('UNINTEGRATED','UNMAPPED'),1])
head(kopa)





## filtrating only those ko by ID_ref mean is > 0.2
sigko <- c()
for (i in levels(sample_ref)){
  sample <- apply(as.matrix(kocov.r.pa.lefe[4:nrow(kocov.r.pa.lefe),which(kocov.r.pa.lefe[1,]==i)]),2,as.numeric)
  sigko <- c(sigko,kocov.r.pa.lefe[4:nrow(kocov.r.pa.lefe),1][which(rowMeans2(sample)>th)])
}  
sigko <- unique(sigko)
sigko <- sigko[!sigko%in%c('UNINTEGRATED','UNMAPPED')]
sigko.lefse <- kocov.r.pa.lefe[kocov.r.pa.lefe$header%in%c(kocov.r.pa.lefe$header[1:3],sigko),]
write.table(sigko,paste(ra,'humann/sigKO.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)
write.table(kocov.r.pa.lefe,paste(ra,'humann/humann2_coverage_ko_lefe.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)
write.table(kocov.r.pa.lefe[1:4,],paste(ra,'humann/metadata_simple.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)


#################### ****###########


### first tests - humann folder ### 

## data for lefse##
tpdata <- as.data.frame(t(pData(ggdata)),stringsAsFactors=FALSE)
tpdata$header <- as.character(rownames(tpdata))
tpdata <- tpdata[,c(ncol(tpdata),1:(ncol(tpdata)-1))]
ko2lfse.r <- read.table(paste(ra,'humann/humann2_pathabundance_relab_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F,skip=2)
colnames(ko2lfse.r) <- colnames(tpdata)
ko2lfse <- rbind.data.frame(tpdata[c('Punto','ID_ref','Fe','Arena','Salinity'),],ko2lfse,stringsAsFactors=F)
write.table(ko2lfse,paste(ra,'humann/humann2_pathabundance_relab_fixed_lefe.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)
write.table(ko2lfse[1:5,],paste(ra,'humann/metadata_simple.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)
## for lefse end ##

### Pathway plot using humann table ##
packages(c('pathview','gage','gageData'))

sample_ref <- factor(kocovsig[1,],levels=c('low','med','high'))

coverage.r <- read.table(paste(ra,'humann/humann2_pathcoverage_contrib_module.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
head(coverage.r)
colnames(coverage.r) <- gsub('_Coverage','',colnames(coverage.r),fixed=T)
tokeep.coverage <- unlist(lapply(as.character(coverage.r[,1]),function(x) length(grep('|',x,fixed=T))==0))
coverage.pa <- coverage.r[tokeep.coverage,]
head(coverage.pa)
tpdata <- as.data.frame(t(pData(ggdata)),stringsAsFactors=FALSE)
tpdata$header <- as.character(rownames(tpdata))
tpdata <- tpdata[,c(ncol(tpdata),1:(ncol(tpdata)-1))]
coverage.pa.lefe <- coverage.pa
colnames(coverage.pa.lefe) <- colnames(tpdata)
coverage.pa.lefe <- rbind.data.frame(tpdata[c('ID_ref','Fe','Arena'),],coverage.pa.lefe,stringsAsFactors=F)
coverage.pa.lefe <- coverage.pa.lefe[-nrow(coverage.pa.lefe),]
## filtrating only those ko by ID_ref mean is > 0.5
th <- 0.5
sigmo <- c()
for (i in levels(sample_ref)){
  sample <- apply(as.matrix(coverage.pa.lefe[4:nrow(coverage.pa.lefe),which(coverage.pa.lefe[1,]==i)]),2,as.numeric)
  sigmo <- c(sigmo,coverage.pa.lefe[4:nrow(coverage.pa.lefe),1][which(rowMeans2(sample)>th)])
}  
sigmo <- unique(sigmo)
sigmo <- sigmo[!sigmo%in%c('UNINTEGRATED','UNMAPPED')]
sigmo.lefse <- coverage.pa.lefe[coverage.pa.lefe$header%in%c(coverage.pa.lefe$header[1:3],sigmo),]
write.table(sigmo,paste(ra,'humann/sigMO.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)
write.table(coverage.pa.lefe,paste(ra,'humann/humann2_coverage_modules_lefse.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)
write.table(coverage.pa.lefe[1:4,],paste(ra,'humann/metadata_simple.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)




kocov.r <- read.table(paste(ra,'humann/humann2_pathcoverage_contrib_Xmp_fixed.txt',sep=''),header=T,sep="\t",quote="",comment.char='*',check.names=F)
#head(kocov.r)
#colnames(kocov.r) <- gsub('_Coverage','',colnames(coverage.r),fixed=T)
tokeep.coverage <- unlist(lapply(as.character(kocov.r[,1]),function(x) length(grep('|',x,fixed=T))==0))
kocov.r.pa <- kocov.r[tokeep.coverage,]
tpdata <- as.data.frame(t(pData(ggdata)),stringsAsFactors=FALSE)
tpdata$header <- as.character(rownames(tpdata))
tpdata <- tpdata[,c(ncol(tpdata),1:(ncol(tpdata)-1))]
kocov.r.pa.lefe <- kocov.r.pa[-c(1,2),]
colnames(kocov.r.pa.lefe) <- colnames(tpdata)
kocov.r.pa.lefe <- rbind.data.frame(tpdata[c('ID_ref','Fe','Arena'),],kocov.r.pa.lefe,stringsAsFactors=F)
kocov.r.pa.lefe <- kocov.r.pa.lefe[-nrow(kocov.r.pa.lefe),]
## filtrating only those ko by ID_ref mean is > 0.2
sigko <- c()
for (i in levels(sample_ref)){
  sample <- apply(as.matrix(kocov.r.pa.lefe[4:nrow(kocov.r.pa.lefe),which(kocov.r.pa.lefe[1,]==i)]),2,as.numeric)
  sigko <- c(sigko,kocov.r.pa.lefe[4:nrow(kocov.r.pa.lefe),1][which(rowMeans2(sample)>th)])
}  
sigko <- unique(sigko)
sigko <- sigko[!sigko%in%c('UNINTEGRATED','UNMAPPED')]
sigko.lefse <- kocov.r.pa.lefe[kocov.r.pa.lefe$header%in%c(kocov.r.pa.lefe$header[1:3],sigko),]
write.table(sigko,paste(ra,'humann/sigKO.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)
write.table(kocov.r.pa.lefe,paste(ra,'humann/humann2_coverage_ko_lefe.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)
write.table(kocov.r.pa.lefe[1:4,],paste(ra,'humann/metadata_simple.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=T)


head(kocov.r.pa.lefe)

dim(kocov.r.pa.lefe)

hist(kocov.r.pa[-c(1,2),])


coverage.p <- coverage.pa[coverage.pa$`2A_Coverage`>0.5,]
coverage.p <- coverage.p[-c(1,2),]
konum <- gsub('ko','',as.character(coverage.p[,1]),fixed=T)
write.table(konum,paste(ra,'humann/KOmaps.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)




colnames(ko2lfse.r) <- colnames(tpdata)

kos <- rep(1,nrow(ko2lfse.r))
names(kos) <- toupper(as.character(ko2lfse.r$header))
head(kos)
kosn <- as.character(ko2lfse.r$header)
komaps <- gsub('ko','map',kosn,fixed=T)
konum <- komaps <- gsub('ko','',kosn,fixed=T)
write.table(konum,paste(ra,'humann/KOmaps.txt',sep=''),quote=F,sep="\t",row.names=F,col.names=F)
ko <- sim.mol.data(mol.type="gene.ko",nmol=5000)

l <- unlist(lapply(kegg.sets.ko,function(x) return(is.element(kosn[1],x))))
table(l)

packages(c('KEGGREST'))

comps <- c()
c <- 0
pb <- txtProgressBar(min = 0, max = length(kosn), style = 3)
for (k in kosn){
  ǵetcomp <- try(comps <- c(comps,keggGet(k)[[1]]$COMPOUND,silent=T))
  if("try-error" %in% class(ǵetcomp)) print(k)
  c <- c+1
  setTxtProgressBar(pb, c)
}

keggGet('KO00312')

pv.out1 <- pathview(gene.data=kosn,cpd.data=unique(names(comps)),species='ko',out.suffix='ko.data',pathway.id='01120',kegg.native=T)

pv.out2 <- pathview(gene.data=names(ko),species='ko',out.suffix='ko.data',pathway.id='04112',kegg.native=T)

node.info(as.character(ko2lfse.r$header))
###



### contributions ##
packages(c('btools'))

diff_genus <- read.table(opt$diff,header=T,sep="\t",quote="")$'Genus'
df_cont <- read.table(opt$contribution,header=T,sep="\t",quote="")
head(df_cont)
contrib <- analyze_contributions(opt$contribution,paste(ra,'metadata_simple.txt',sep=''))
#max(df_cont$ContributionPercentOfSample)
length(levels(df_cont$Genus))
k <- df_cont[grep(as.character(diff_genus[4]),df_cont$Genus),]
k <- k[with(k,order(-ContributionPercentOfAllSamples)),]

#dendjw <- stringdistmatrix(strsplit(as.character(diff_genus[4]),'_')[[1]][1],df_cont$Genus,method='jw',useNames='strings')
#dendjw[1,c(151,162)]
#which(dendjw==max(dendjw))

#%>%as.dist %>% hclust(method='average') %>% as.dendrogram


ggplot(k,aes(x=Sample,y=ContributionPercentOfAllSamples,fill=Genus))+geom_boxplot()



