#!/usr/bin/env Rscript
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
# Rscript beta_div.m.R -h
# Rscript beta_div.m.R -i ~/16Srlib_test/results/dataF.rds -o ~/16Srlib_test/results/ -e [options] 
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
packages(c("optparse"))

## Options ##

#p <- '/home/torres/ikmb_storage/Mangrove/16Sfa/08_2017_results/2017/'
p <- '/Users/guillermotorres/Documents/Proyectos/UAN_Manglares/29_08_ResultsMgv16s_v2/'

option_list <- list(
  make_option(c("-i","--data"),action="store",type="character",default=paste(p,'dataFcp_l_0.1.rds',sep=''),
              help="Path to input rds file"),
  make_option(c("-o","--out"),action="store",type="character",default=paste(p,'beta/',sep=''),
              help="Path to output directory [default %default]"),
  make_option(c("-e","--exploratory"),action="store_true",default="NA",
              help="Perform exploratory analysis"),
  make_option(c("-m","--model"),action="store_true",default="NA",
              help="Build constrained model based on AIC selection criterion"),
  make_option(c("-am","--assess_model"),action="store",type="character",default=NA,#NULL,
              help="Model's terms assessed by permutation tests; for new model: vars,separated,by,comma"),
  make_option(c("-C","--constraints"),action="store",type="character",default=NULL,
              help="Set of constraints used by -coa: vars,separated,by,comma"),
  make_option(c("-coa","--constrained_analysis"),action="store_true",default="NA",
              help="Perform constrained ordination analysis using -C constraints"),
  make_option(c("-a","--factors"),action="store",type="character",default=NULL,
              help="Set of factors used by -b: vars,separated,by,comma"),
  make_option(c("-b","--beta"),action="store_true",default='NA',
              help="Perform beta diversty between -a variable classes"),
  make_option(c("-f","--filter"),action="store",type="double",default=0.3,
              help="Percentile as threshold of low abundant features ")
)
parser <- OptionParser(usage = "%prog -i path/to/infile -o path/to/outdir [options]",option_list=option_list)
opt <- parse_args(parser)
#parse_args(parser,positional_arguments=1) 
if (is.na(opt$data)){stop(sprintf("There is not file specified"))}

##
if(dir.exists(opt$out)){message('Out-folder already exist, files will be overwritten')
}else dir.create(opt$out,showWarnings=F)
#### Preparing the input data ####
packages(c("metagenomeSeq","reshape2","vegan","ggplot2"))
data.raw <- readRDS(opt$data)
pData(data.raw) <- pData(data.raw)[, colSums(is.na(pData(data.raw))) != nrow(pData(data.raw))] # rm columns with all values nas
quant <- pData(data.raw)[, !sapply(pData(data.raw), is.factor)] ## quantitative explanatory variables 
cate <- pData(data.raw)[, sapply(pData(data.raw), is.factor)] ## categoric explanatory variables

#if (!"ID"%in%colnames(pData(data))) pData(data)$ID <- pData(dfp)[,1]
## taxonomical level ##
tl='Genus'
data <- aggTax(data.raw,lvl=tl)

# Data transformations 
dfr <- t(MRcounts(data)) # raw data
dfn <-t(MRcounts(data,norm=T)) # cumsum normalization
dfnl <-t(MRcounts(data,norm=T,log=T)) # cumsum normalization plus log transform
dfhel <-decostand(t(MRcounts(data)),'hellinger') ## Hellinger normalization
dflhel <-decostand(t(MRcounts(data,log=T)),'hellinger') ## log data Hellinger normalization
dfchor <-decostand(t(MRcounts(data)),'norm') ## log data Hellinger normalization

####### end ######

### Forward selection of explanatory variables ##
## search for parsimony, and possible strong linear dependencies (correlations) among the explanatory variables in the RDA model

## Linear dependencies can be explored by computing the variables’ variance inflation factors (VIF), which measure the proportion by which the variance of a regression coefficient is inflated in the presence of other explanatory variables.
sink(file=paste(opt$out,"forwardSelection_R2_RDA.txt",sep='')) 
forward.sel <- ordiR2step(rda(dfhel~1,data=quant),scope=formula(rda(dfhel~.,data=quant)),R2scope=F,direction='forward',pstep=1000,R2permutations=1000)
sink(NULL)
sigVars <- attr(forward.sel$terms,"term.labels")

### RDA calculation ####
df <-  'dfhel'
df.rda <- rda(formula(paste(df,'~',paste(sigVars,collapse='+'),sep='')),data=quant)
# R2 measure the unibiased amount of explained variation. R-squared is a statistical measure of how close the data are to the fitted regression line. R2 is the percentage of the response variable variation that is explained by a linear model. R-squared = Explained variation / Total variation
# Ezekiel’s adjustment (Ezekiel 1930). As a rule of thumb, this adjustment may be overly conservative when m>n/2. Where n is the number of objects (samples) and m is the number of explanatory variables.
R2 <- RsquareAdj(df.rda)
rdasummary <- summary(df.rda)

## writing information of RDA ##
sink(file=paste(opt$out,"summary_RDA.txt",sep='')) 
print('Linear model coefficients from explanatory variables:')
coef(df.rda)
print('R2 is the percentage of the response variable variation that is explained by a linear model: Calculated and Ezekiel’s adjustment :')
R2
rdasummary
sink(NULL)

## Variance represented #
names(rdasummary)
# Overall variance - corrected by R2
constvar <- sum(rdasummary$concont$importance[2,]*R2$adj.r.squared)
#  The constrained fraction is the amount of variance of the Y matrix explained by the explanatory variables. Expressed as a proportion, it is equivalent to an R2 in multiple regression
unconstvar <- 1-constvar
#* The canonical (RDAx) eigenvalues measure amounts of variance explained by the RDA model, whereas the residual (PCx) eigenvalues measure amounts of variance represented by the residual axes, but not explained by the RDA model.
#
RDA1cont <- R2$adj.r.squared*rdasummary$concont$importance[2,1]
RDA2cont <- R2$adj.r.squared*rdasummary$concont$importance[2,2]
RDA3cont <- R2$adj.r.squared*rdasummary$concont$importance[2,3]

# Explanatory Variables explained variance:

expVar_explained <- rdasummary$concont$importance[2,]*R2$adj.r.squared
names(expVar_explained) <- rownames(rdasummary$biplot)

### Permutation test to evaluate our RDA results.
# Global test first
anova.cca(df.rda,step=1000)
# test of all canonical axes
anova.cca(df.rda,by='axis',step=1000)
anova.cca(df.rda,by='terms',step=1000)
anova.cca(df.rda,by='margin',step=1000)

sink(file=paste(opt$out,"RDA_importantInfo.txt",sep=''))
paste('Overall variance corrected by R2; partitioned into constrained:',sep='')
constvar
paste('and unconstrained:')
unconstvar
paste('Contribution of the explanatory variables: ',sep='')
expVar_explained
paste('Permutation test. Validating our general model: ',sep='')
anova.cca(df.rda,step=1000)
paste('Permutation test. Validating each canonical axes: ',sep='')
anova.cca(df.rda,by='axis',step=1000)
paste('Permutation test. Validating each explanatory variables: ',sep='')
anova.cca(df.rda,by='terms',step=1000)
sink(NULL)

## wa scores
#scaling both 
# scaling 1: distance biplot: 1) Projecting an object (sample) at right angle on a response variable (species) or a quantitative explanatory variable approximates the position of the object along that variable. 2) The angles between response (species) and explana- tory variables in the biplot reflect their correlations 3) The relationship between the centroid of a qualitative explanatory variable and a response variable (species) is found by projecting the centroid at right angle on the species variable 4) Distances among centroids, and between centroids and individual objects, approximate their Euclidean distances.
# scaling 2: correlation biplot: 1) Projecting an object (sample) at right angle on a response variable (species) or a quantitative explanatory variable approximates the value of the object along that variable. 2) The angles in the biplot between response and explanatory variables, and between response variables themselves or explana- tory variables themselves, reflect their correlations 3) The relationship between the centroid of a qualitative explanatory variable and a response vari- able (species) is found by projecting the centroid at right angle on the species variable (as for individual objects) 4) Distances among centroids, and between centroids and individual objects, do not approximate their Euclidean distances.

sp.scores <- scores(df.rda,choices=1:2,scaling=1,display='sp')
site.socres <- scores(df.rda,choices=1:2,scaling=1,display='site')
newnames <- data.frame(oldnames=rownames(MRcounts(data)),newnames=make.cepnames(rownames(MRcounts(data))))


packages(c('Rlof',"RColorBrewer","RAM",'directlabels'))

## colors of the sites ###
extravar <- 'ID_ref'
pData(data)[[extravar]] <- factor(pData(data)[[extravar]],levels=c('low','med','high'))
index <- pData(data)
if (length(levels(index[[extravar]]))>9){ l <- length(levels(index[[extravar]]))
}else l <- 9
sampleColors= colorRampPalette(brewer.pal(8, "Set1"))(l)[c(1:length(levels(index[[extravar]])))]
sample_colors_dframe = data.frame(samples=levels(index[[extravar]]), colors=sampleColors)
sample_colors = unlist(lapply(sort(index[[extravar]]), function(x) return(as.character(sample_colors_dframe[x,2]))))
names_data <- data.frame(rowids=colnames(data),ref=index[match(as.character(index[[1]]),colnames(data))][[extravar]])
names_data <- names_data[with(names_data, order(ref)),]
names_data$colors <- sample_colors
## end colors ###

## identify external OTUs using outlier approach
## outliers using Local Outlier Factor (LOF) approach##
l <- round(lof(sp.scores,nrow(sp.scores)/2),1)
x <- 1.3
outliers <- which(l>x)
cols <- rep("grey",length(l))
cols[outliers]='red'
### end outliers ## 

# Plot scaling 1:
pdf(paste(opt$out,"RDA_scaling1_2more_names.pdf",sep=''),width=9, height=8,bg=T,useDingbats=F)
fig <- ordiplot(df.rda,type="none",display=c('sites','sp'),scaling=1)
text(df.rda,display='cn',col='blue',scaling=1,axis.bp=T)
points(fig,"sp",pch=3,col=cols,cex=0.8)
points(site.socres,pch=16,cex=0.9,col=names_data$colors[match(rownames(site.socres),names_data$rowids)])
ordihull(df.rda,groups=names_data$ref[match(rownames(site.socres),names_data$rowids)],display='sites',scaling=1,draw='polygon',label=T,col=sample_colors_dframe$colors)
#text(sp.scores[outliers,],labels=rownames(sp.scores)[outliers],pos=3,cex=0.7)
#textplot(sp.scores[outliers,1],sp.scores[outliers,2],words=newnames$newnames[outliers],cex=0.8,new=F)
orditorp(df.rda,display='sp',select=c(l>x),labels=as.character(newnames$newnames[outliers]),scaling=1,cex=0.8,pch=3,pcol='red',air=0.1,pos=3)
#orditorp(df.rda,display='sp',select=c(l>x),labels=rownames(sp.scores)[outliers],scaling=1,cex=0.8,pch=3,pcol='red',air=0.1,pos=3)
dev.off()

write.table(newnames[outliers,],paste(opt$out,"RDAplot_names.txt",sep=''),quote=F,sep="\t",row.names=F)


plot(df.rda,scaling=1,main='Triplot RDA scaling 1 - wa scores',display=c('sp','wa','cn'),xlab=paste("RDA1: ",round(RDA1cont*100),'%',sep=''),
     ylab=paste("RDA2: ",round(RDA2cont*100),'%',sep=''))
#text(df.rda.scores, labels=rownames(df.rda.scores), pos=3)

arrows(0,0,df.rda.scores[,1],df.rda.scores[,2],length=0,lty=1,col='red')

plot(df.rda,scaling=1,main='Triplot RDA scaling 1 - wa scores',display=c('sp','wa','cn'),choices=c(1,3))


### Scaling 2
## identify external OTUs using outlier approach
## outliers using Local Outlier Factor (LOF) approach##
l <- round(lof(sp.scores1,nrow(sp.scores)/2),1)
x <- 1.5
outliers <- which(l>x)
cols <- rep("grey",length(l))
cols[outliers]='red'
### end outliers ## 

sp.scores2 <- scores(df.rda,choices=1:2,scaling=2,display='sp')
site.socres2 <- scores(df.rda,choices=1:2,scaling=2,display='site')

pdf(paste(opt$out,"RDA_scaling2.pdf",sep=''),width=9, height=9,bg=T,useDingbats=F)
#fig <- ordiplot(df.rda,type="none",display=c('sp','sites'),scaling=2,xlim=c(-0.5,0.5),ylim=c(-0.5,0.7))
fig <- ordiplot(df.rda,type="none",display=c('sp'),scaling=2)
text(df.rda,display='cn',col='blue',scaling=2,axis.bp=T)
points(fig,"sp",pch=3,col=cols,cex=0.8)
points(site.socres2,pch=16,cex=0.9,col=names_data$colors[match(rownames(site.socres),names_data$rowids)])
ordihull(df.rda,groups=names_data$ref[match(rownames(site.socres),names_data$rowids)],display='sites',scaling=2,draw='polygon',label=T,col=sample_colors_dframe$colors)
orditorp(df.rda,display='sp',select=c(l>x),labels=as.character(newnames$newnames[outliers]),scaling=2,cex=0.8,pch=3,pcol='red',air=0.1,pos=3)
#orditorp(df.rda,display='sp',select=c(l>x),labels=rownames(sp.scores)[outliers],scaling=1,cex=0.8,pch=3,pcol='red',air=0.1,pos=3)
dev.off()
write.table(newnames[outliers,],paste(opt$out,"RDA_s2_plot_names.txt",sep=''),quote=F,sep="\t",row.names=F)


plot(df.rda,main='Triplot RDA scaling 2 - wa scores',display=c('sp'),scaling=2)
text(df.rda,display='cn',col='blue',scaling=2,axis.bp=T)
orditorp(df.rda,display='sp',select=c(l>x),labels=as.character(newnames$newnames[outliers]),scaling=2,cex=0.8,pch=3,pcol='red',air=0.1,pos=3)



## lc scores

plot(df.rda,scaling=1,main='Triplot RDA scaling 1 - lc scores',display=c('sp','lc','cn'))
df.rda.scores <- scores(df.rda,choices=1:2,scaling=1,display='sp')
arrows(0,0,df.rda.scores[,1],df.rda.scores[,2],length=0,lty=1,col='red')

plot(df.rda,main='Triplot RDA scaling 2 - wa scores',display=c('sp','lc','cn'))
df.rda.scores2 <- scores(df.rda,choices=1:2,display='sp')
arrows(0,0,df.rda.scores2[,1],df.rda.scores2[,2],length=0,lty=1,col='red')


## Apply kaiser-gutmar criterion to residual axes

df.rda$CA$eig[df.rda$CA$eig>mean(df.rda$CA$eig)]


### RDA as a tool for multivariate anova
## salinity
dfhel.MHV <- betadisper(dist(dfhel),pData(data)[,2])
anova(dfhel.MHV)
permutest(dfhel.MHV)

## Tukey's Honest Significant Differences
dfhel.HSD <- TukeyHSD(dfhel.MHV)
plot(dfhel.HSD)
plot(dfhel.MHV, ellipse = TRUE, hull = FALSE, conf = 0.90)

## soil structure 

dfhel.ss.MHV <- betadisper(vegdist(dfn),pData(data)$Textura)
anova(dfhel.ss.MHV)
permutest(dfhel.ss.MHV)

## Tukey's Honest Significant Differences
dfhel.ss.HSD <- TukeyHSD(dfhel.ss.MHV)
plot(dfhel.ss.HSD)
plot(dfhel.ss.MHV, ellipse = TRUE, hull = FALSE, conf = 0.90)


#interaction.rda <- rda(dfhel,,)


