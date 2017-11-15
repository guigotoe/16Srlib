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
pa <- '/home/torres/ikmb_storage/Mangrove/16Sfa/2016/mothur/'                   # MANGROVE VERSION
ra <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2016/'          # MANGROVE VERSION

option_list <- list(
  make_option(c("-i","--input"),action="store",type="character",default=paste(pa,'greengenes_result/16s_2016_salinity.biom',sep=''),#'16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared'
              help="Path to input biom input file"),
  make_option(c("-m","--metadata"),action="store",type="character",default=paste(pa,'metadata.txt',sep=''),
              help="Path to input metadata file"),
  make_option(c("-f","--fasta"),action="store",type="character",default=paste(pa,'OTUs.fa',sep=''),
              help="Path to OTUs fasta file"),
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
picrust_outpath <- paste(opt$out,'picrust/',sep='')
human_outpath <- '/home/sukmb347/sukmb347/Mangrove/16Sfa/08_2017_results/2016/human/'# **server version 

if(dir.exists(picrust_outpath)){message('picrustOut folder already exist, files will be overwritten')
}else dir.create(picrust_outpath,showWarnings=F)

###### end ######
packages(c('biomformat','phyloseq'))
ggdata <- mothur.biom(opt$input,opt$metadata)
p.f <- cumNormStatFast(ggdata) # Calculates the percentile for which to sum counts up to and scale by.
ggdata <- cumNorm(ggdata,p=p.f)  # Calculates each column's quantile and calculates the sum up to and including p quantile
nf.f <- normFactors(ggdata)
head(fData(ggdata))
##* This biom export option is not working correctly *##
ggbiom <- MRexperiment2biom(ggdata)
ggbiom <- make_biom(round(MRcounts(ggdata,norm=T)),sample_metadata=pData(ggdata),observation_metadata=fData(ggdata))#,sample_metadata=pData(ggdata),observation_metadata=fData(ggdata),id='Mangrove16s2016')

newbiom <- paste('CumSum_',rev(strsplit(opt$input,'/',fixed=T)[[1]])[1],sep='')
write_biom(biom(ggbiom),paste(picrust_outpath,newbiom,sep=''))

showMethods(ggbiom)

read_biom(biom(ggbiom))
read_biom(ggbiom)

write_biom(biom(MRexperiment2biom(ggdata)),paste(picrust_outpath,newbiom,sep=''))
### end ###

## * PICRUST * ##
CpNnormalization <- paste('python ',picrust_path,'normalize_by_copy_number.py -i ',picrust_outpath,newbiom,' -o ',
                          picrust_outpath,'normalized_otus.biom -c /home/torres/Bin/picrust-1.1.2/picrust/data/16S_13_5_precalculated.tab.gz',sep='')
system(CpNnormalization)
## KO prediction : calculates CI and Normalized the predicted functional abundances by dividing each abundance by the sum of functional abundances in the sample. Total sum of abundances for each sample will equal 1.
KOprediction <- paste('python ',picrust_path,'predict_metagenomes.py -f -i ',picrust_outpath,'normalized_otus.biom -o ',
                      picrust_outpath,'ko_predicted_metagenomes.txt --with_confidence --normalize_by_function --suppress_subset_loading',sep='')
system(KOprediction)
KOprediction_biom <- paste('python ',picrust_path,'predict_metagenomes.py -i ',picrust_outpath,'normalized_otus.biom -o ',
                      picrust_outpath,'ko_biom_predicted_metagenomes.biom --with_confidence --normalize_by_function --suppress_subset_loading',sep='')
system(KOprediction_biom)

## COG prediction
COGprediction <- paste('python ',picrust_path,'predict_metagenomes.py -f -i ',picrust_outpath,'normalized_otus.biom -o ',
                       picrust_outpath,'cog_predicted_metagenomes.txt -t cog --with_confidence --normalize_by_function --suppress_subset_loading',sep='')
system(COGprediction)
## rfam prediction
Rfamprediction <- paste('python ',picrust_path,'predict_metagenomes.py -f -i ',picrust_outpath,'normalized_otus.biom -o ',
                        picrust_outpath,'rfam_predicted_metagenomes.txt -t rfam --with_confidence --normalize_by_function --suppress_subset_loading',sep='')
system(Rfamprediction)
# Output is a tab-delimited column indicating OTU contribution to each function.
ko_contributions <- paste('python ',picrust_path,'metagenome_contributions.py -i ',picrust_outpath,'normalized_otus.biom -o ',
                        picrust_outpath,'ko_contributions.txt -t ko --suppress_subset_loading',sep='')
system(ko_contributions)
COG_contributions <- paste('python ',picrust_path,'metagenome_contributions.py -i ',picrust_outpath,'normalized_otus.biom -o ',
                          picrust_outpath,'cog_contributions.txt -t cog --suppress_subset_loading',sep='')
system(COG_contributions)
Rfam_contributions <- paste('python ',picrust_path,'metagenome_contributions.py -i ',picrust_outpath,'normalized_otus.biom -o ',
                          picrust_outpath,'rfam_contributions.txt -t rfam --suppress_subset_loading',sep='')
system(Rfam_contributions)

## collapse 

ko_colapse_l1 <- paste('python ',picrust_path,'categorize_by_function.py -i ',picrust_outpath,'ko_biom_predicted_metagenomes.biom -o ',
                       picrust_outpath,'ko_predicted_metagenomes.L1.txt -f -l 1 -c KEGG_Pathways',sep='')
system(ko_colapse_l1)
ko_colapse_l1.b <- paste('python ',picrust_path,'categorize_by_function.py -i ',picrust_outpath,'ko_biom_predicted_metagenomes.biom -o ',
                       picrust_outpath,'ko_predicted_metagenomes.L1.biom -l 1 -c KEGG_Pathways',sep='')
system(ko_colapse_l1.b)

ko_colapse_l2 <- paste('python ',picrust_path,'categorize_by_function.py -i ',picrust_outpath,'ko_biom_predicted_metagenomes.biom -o ',
                       picrust_outpath,'ko_predicted_metagenomes.L2.txt -f -l 2 -c KEGG_Pathways',sep='')
system(ko_colapse_l2)

ko_colapse_l3 <- paste('python ',picrust_path,'categorize_by_function.py -i ',picrust_outpath,'ko_biom_predicted_metagenomes.biom -o ',
                       picrust_outpath,'ko_predicted_metagenomes.L3.txt -f -l 3 -c KEGG_Pathways',sep='')
system(ko_colapse_l3)
##

#ko_colapse_tax <- paste('python ',picrust_path,'categorize_by_function.py -i ',picrust_outpath,'ko_biom_predicted_metagenomes.biom -o ',
#                       picrust_outpath,'ko_genus_table.txt -f -l 6 -c taxonomy',sep='')
#system(ko_colapse_tax)

#


# HUMAnN2 (the HMP Unified Metabolic Analysis Network) is a method for efficiently and accurately determining the 
# presence, absence, and abundance of metabolic pathways in a microbial community from metagenomic or metatranscriptomic sequencing data
# ** Run HUMAnN2 in server.. everithing is already installed

## Constructing the exec-file to use on the server.

header <- paste('#!/bin/bash
#SBATCH --job-name=micro_l
#SBATCH --output=microl.out
#SBATCH --error=microl.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120000
#SBATCH --partition=msb
#SBATCH --qos=msb
#SBATCH --mail-user=g.torres@ikmb.uni-kiel.de
#SBATCH --mail-type=ALL
export OMP_NUM_THREADS=1 \ncd ',human_outpath,
'\nmodule load IKMB
module load Python/2.7.14
module load Bowtie2/2.3.2
module load Diamond/0.7.9
module load Minpath/1.2
module load Metaphlan/2.0 
module load Humann2/0.6.1',sep='')
write(header,paste(opt$out,"human/human.exec",sep=''))







