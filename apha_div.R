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
# Rscript apha_div.R /path/dataF.rds
#
#* requirements *#

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "toolbox.R")
source("/home/torres/Documents/Projects/Metagenome/r_scripts/16Srlib/toolbox.R")
source(other.name)

packages(c("metagenomeSeq","vegan","ggplot2"))

###### end ######

#* input *

f <- '/home/torres/Documents/Projects/Metagenome/results/plotsMothur/09.2016/dataF.rds' # commandArgs()[6] #
df <- readRDS(f)

## diversity##

shannon <- vegan::diversity(t(MRcounts(df,norm=T)))
shannon <- as.data.frame(shannon)
shannon$id <- rownames(shannon)
index <- merge(pData(df),shannon,by="id")
assocvar <- "Gender"
xvar <- "Age"


ggplot(index,aes(x=xvar,y=shannon,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
  scale_colour_hue(name=Gender,breaks=levels(index$assocvar),labels=clevels(index$assocvar),l=50)+
  stat_smooth(method="lm")+facet_grid(.~Gender,margins=TRUE)

## for Age ##
#index <- merge(pData(df),shannon,by="id")
#index_sum <- summarySE(index[,c("age.mean","Gender","shannon")],measurevar="shannon",groupvars=c("age.mean","Gender"),na.rm=FALSE)

#ggplot(index[!is.na(index$Gender),],aes(x=Age,y=shannon,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
#  xlab("Individuals age")+scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+
#  stat_smooth(method="lm")+facet_grid(.~Gender,margins=TRUE)
# 
#ggplot(index_sum,aes(x=age.mean,y=shannon,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
#  xlab("Average group's age")+ scale_colour_hue(name="Gender",l=50)+stat_smooth(method="lm")#+







D1 <- summarySE(df[,c("age.mean","Gender","D1")],measurevar="D1",groupvars=c("age.mean","Gender"),na.rm=FALSE)
ggplot(df[!is.na(df$Gender),],aes(x=Age,y=D1,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
  xlab("Individuals age")+scale_colour_hue(name="Gender",breaks=c("female", "male"),labels=c("Female", "Male"),l=50)+
  stat_smooth(method="lm")+facet_grid(.~Gender,margins=TRUE)
if(!is.null(savef)) ggsave(paste(savef,title,'_shannonD1.pdf',sep=""),width=12, height=8)
ggplot(D1,aes(x=age.mean,y=D1,col=Gender))+geom_point()+ylab("Alpha diversity (Shannon entropy)")+
  xlab("Average group's age")+
  scale_colour_hue(name="Gender",l=50)+stat_smooth(method="lm")#+

