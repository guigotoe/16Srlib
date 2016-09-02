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
####################################################
# Prepare and filters the data from mothur.
# How to use:
# Rscript data_prep.R /path/counts_file /path/taxonomy_file /path/metadata.txt min_percentage_OTU_presence[0-1] /path/for/out_file/
# Rscript data_prep.R /path/16s.an.shared /path/16s.an.cons.taxonomy /path/metadata.txt 0.05 /path/for/out_file/
#
#* requirements *#

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "toolbox.R")
source(other.name)
#source("/home/torres/Documents/Projects/Metagenome/bin/rscripts/16Srlib/toolbox.R")
#print(paste("Sourcing",other.name,"from",script.name))
packages(c("metagenomeSeq"))

###### end ######

#* input *
c <- commandArgs()[6] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/16s.an.shared' # commandArgs()[6] #
t <- commandArgs()[7] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/16s.an.cons.taxonomy' # commandArgs()[7] #
m <- commandArgs()[8] #'/home/torres/ikmb_storage/Metagenome/16s/03.2016/metadata.txt' # commandArgs()[8] #
th <- commandArgs()[9] #0.05 #commandArgs()[9] # percentage threshold of OTU's presence across the samples.
o <- commandArgs()[10] #'/home/torres/Documents/Projects/Metagenome/results/plotsMothur/09.2016/' #
d <- 10 # depth count threshold - by default.
  
message("Preparing the files...")
metadata <- mothur.metadata(read.table(m,header=T,sep="\t",row.names=1,blank.lines.skip=TRUE,na.strings=c("","NA")))
taxonomy <- mothur.taxonomy(read.table(t,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
counts <- mothur.counts(read.table(c,header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA")))
ord = match(colnames(counts),rownames(metadata))
metadata = metadata[ord,]
data <- newMRexperiment(counts,phenoData=AnnotatedDataFrame(metadata),featureData=AnnotatedDataFrame(taxonomy))
message("Filtering...")
data.f <- filterData(data,present=round(th*NROW(pData(data))),depth=d)
retained.info <- sum(fData(data.f)$Size)/sum(fData(data)$Size)
message(paste("OTUs with less than ",d," counts and their presence in less than ",100*th,"% of the samples, were removed\n",
              "From ",dim(MRcounts(data))[1]," OTUs to ",dim(MRcounts(data.f))[1],"\n",
              "Information retained: ",round(100*retained.info,2),"%; lost: ",round((100-100*retained.info),2),"%",sep=""))
## Normalizing
message("Normalizing...")
p.f <- cumNormStatFast(data.f) # Calculates the percentile for which to sum counts up to and scale by.
data.f <- cumNorm(data.f,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(data.f)
# saving files
message("Exporting files...")
saveRDS(data.f,file=paste(o,'dataF.rds',sep=''))
message("dataF.rds -> Counts filtered and normalized - Scaling Factor Normalization (Paulson et. al 2013)\n",
        "Data was successfully prepared!")


