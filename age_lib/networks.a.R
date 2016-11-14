####################################################
# By Guillermo Torres PhD.c                        #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #            
####################################################
# Last update: 05.08.2016
# Created: 04 August 2016
#
# This is written as part of 16S - Aging analysis, but could be
# splitted to serve different purposes.
####################################################

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
#toolbox <- "/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib/toolbox.R"
source(toolbox)

packages(c("igraph","network","sna","ergm","tsna","ndtv","genefilter","visNetwork","networkDynamic"))

## Options ##
p <- '/home/torres/Documents/Projects/Metagenome/MothurResults/03.2016/'
#p <- '/Users/guillermotorres/Documents/Proyectos/Doctorado/16Srlib_test/'


fullnet <- readRDS(paste(p,'fullnet.mb.rds',sep=''))

mnetg <- fullnet

names(mnetg)
names(mnetw)
V(b)$phylla <- as.vector(unlist(unname(taxclass.mc["phylla",])))
V(b)$phycol <- brewer.pal(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["phylla",])))))), "Set1")
V(b)$genus <-as.vector(unlist(unname(taxclass.mc["genus",])))
V(b)$gencol <- colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["genus",])))))))
V(b)$family <- as.vector(unlist(unname(taxclass.mc["family",])))
V(b)$famcol <- colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["family",])))))))
V(b)$type <- "OTU"


phylla_colors <- levels(as.factor(V(mnetg$full$igraph)$phycol))
genus_colors <- levels(as.factor(V(mnetg$full$igraph)$gencol))
family_colors <- levels(as.factor(V(mnetg$full$igraph)$famcol))

par(mfrow=c(2,2))
par(mar = c(1, 1, 2, 15))
 
plot(mnetg[["full"]][["igraph"]], layout=mnetg[["full"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["full"]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(mnetg[["full"]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main="All indiduals - Phylla")
legend(x=1.1, y=1.3, levels(as.factor(V(mnetg[["full"]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(mnetg[["full"]][["igraph"]], layout=mnetg[["full"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["full"]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(mnetg[["full"]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main="All indiduals - Famillies")
legend(x=1.1, y=1.3, levels(as.factor(V(mnetg[["full"]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
plot(mnetg[["full"]][["igraph"]], layout=mnetg[["full"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["full"]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(mnetg[["full"]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main="All indiduals - Genders")
legend(x=1.1, y=1.3, levels(as.factor(V(mnetg[["full"]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=2)
###

plot(mnetg[["30"]][["igraph"]], layout=mnetg[["30"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["30"]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(mnetg[["30"]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["30"]))
legend(x=-1, y=-1.1, levels(as.factor(V(mnetg[["30"]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=4)

plot(mnetg[["50"]][["igraph"]], layout=mnetg[["50"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["50"]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(mnetg[["50"]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["50"]))
legend(x=-1, y=-1.1, levels(as.factor(V(mnetg[["50"]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=4)
plot(mnetg[["70"]][["igraph"]], layout=mnetg[["70"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["70"]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(mnetg[["70"]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["70"]))
legend(x=-1, y=-1.1, levels(as.factor(V(mnetg[["70"]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=4)
plot(mnetg[["85"]][["igraph"]], layout=mnetg[["85"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["85"]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(mnetg[["85"]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["85"]))
legend(x=-1, y=-1.1, levels(as.factor(V(mnetg[["85"]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=4)


## family

plot(mnetg[["30"]][["igraph"]], layout=mnetg[["30"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["30"]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(mnetg[["30"]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["30"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["30"]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=3)
plot(mnetg[["50"]][["igraph"]], layout=mnetg[["50"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["50"]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(mnetg[["50"]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["50"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["50"]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=3)
plot(mnetg[["70"]][["igraph"]], layout=mnetg[["70"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["70"]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(mnetg[["70"]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["70"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["70"]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=3)
plot(mnetg[["85"]][["igraph"]], layout=mnetg[["85"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["85"]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(mnetg[["85"]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["85"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["85"]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=3)

## genus

plot(mnetg[["30"]][["igraph"]], layout=mnetg[["30"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["30"]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(mnetg[["30"]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["30"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["30"]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=4)
plot(mnetg[["50"]][["igraph"]], layout=mnetg[["50"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["50"]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(mnetg[["50"]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["50"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["50"]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=4)
plot(mnetg[["70"]][["igraph"]], layout=mnetg[["70"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["70"]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(mnetg[["70"]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["70"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["70"]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=4)
plot(mnetg[["85"]][["igraph"]], layout=mnetg[["85"]][["coord"]],rescale=T, 
     vertex.size=mnetg[["85"]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(mnetg[["85"]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(mnetg["85"]))
legend(x=-2, y=-1.1, levels(as.factor(V(mnetg[["85"]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=4)


dim(gtensor.m$"30")
dim(gtensor.m$"50")
dim(gtensor.m$"70")
dim(gtensor.m$"85")

######################################################################
mnetw <- readRDS(paste(filesp,'fnet_wl20_40win.mb.rds',sep=''))

netlist <- c()
for (i in names(mnetw)){
  netlist <- c(netlist,mnetw[[as.character(i)]]$se.mb$refit)
}
netlist<-lapply(netlist,as.network.matrix,matrix.type='adjacency')
mnet<-networkDynamic(network.list=netlist)
p.names <- V(mnetw[[as.character(i)]]$igraph)$phylla
f.names <- V(mnetw[[as.character(i)]]$igraph)$family
g.names <- V(mnetw[[as.character(i)]]$igraph)$genus
network.vertex.names(mnet)<-g.names
mnet%v%'family' <- f.names
mnet%v%'phyllum' <- p.names

mnet%n%'net.obs.period'<-list(
  observations=list(c(0,length(names(mnetw)))),
  mode="discrete", 
  time.increment=1,
  time.unit="book volume")

render.animation(mnet)
compute.animation(mnet,
                  animation.mode='MDSJ',
                  default.dist=2,
                  verbose=FALSE)
render.d3movie(mnet,
               render.par=list(tween.frames=20),
               vertex.cex=0.8,label.cex=0.8,label.col='gray',
               displaylabels=F,
               # make shape relate to school year
               #vertex.sides=mnet%v%'schoolyear'-1983,
               # color by gender
               #vertex.col=ifelse(harry_potter_support%v%'gender'==1,'blue','green'),
               edge.col="darkgray",
               output.mode = 'htmlWidget'
)
m <- mnet
par(mfrow=c(1,1))
par(mar = c(3, 3, 2, 4))

proximity.timeline(mnet,default.dist=6,
                   mode='sammon',labels.at=124,vertex.cex=4)
plot(tEdgeFormation(mnet))
plot(tSnaStats(mnet,'gtrans'))
plot(tErgmStats(mnet,'meandeg'),main='Mean degree of Men nets')

timeline(mnet,slice.par=list(start=0,end=25,interval=1,aggregate.dur=1,rule='latest'),
         plot.vertex.spells=FALSE)

timePrism(mnet,at=c(1,10,20),displaylabels=T,planes=T,label.cex=0.5)
timeline(mnet)

######################################### old bases ##################

## Hubs and authorities

hs <- hub_score(tensor.f.mb[[46]][["igraph"]],weights=NA)$vector
as <- authority_score(tensor.f.mb[[46]][["igraph"]],weights=NA)$vector

par(mfrow=c(1,2))
plot(tensor.f.mb[[46]][["igraph"]],vertex.size=hs*50,main="Hubs")
plot(tensor.f.mb[[46]][["igraph"]],vertex.size=as*30,main="Authorities")

## community detection

ceb <- cluster_edge_betweenness(tensor.f.mb[[46]][["igraph"]])
dendPlot(ceb,mode="hclust")
plot(ceb,tensor.f.mb[[46]][["igraph"]])


cfg <- cluster_fast_greedy(as.undirected(tensor.f.mb[[46]][["igraph"]]))
plot(cfg, as.undirected(tensor.f.mb[[46]][["igraph"]]))

colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
kc <- coreness(tensor.f.mb[[46]][["igraph"]], mode="all")
plot(tensor.f.mb[[46]][["igraph"]], vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])

###
par(mfrow=c(1,1))
par(mar = c(5, 5, 1, 1))


par(mfrow=c(1,1))
par(mar = c(0, 0, 1.3, 0))

phylla_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["phylla",])))))))
genus_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["genus",])))))))
family_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(levels(as.factor(as.vector(unlist(unname(taxclass.mc["family",])))))))


V(tensor.f.mb[[46]][["igraph"]])$colphy <- tax_colors_dframe[V(tensor.f.mb[[46]][["igraph"]])$phylla]

shape <- c("square","circle","rectangle","crectangle","vrectangle")

plot(tensor.f.mb[[46]][["igraph"]], layout=tensor.f.mb[[46]][["coord"]],rescale=T, 
     vertex.size=tensor.f.mb[[46]][["vsize"]], 
     vertex.color=phylla_colors[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     #vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     edge.curved=.15,
     edge.width=2,
     main=names(tensor.f.mb[46]))
legend(x=-1, y=-1.1, levels(as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)), pch=21,
       col="#777777", pt.bg=phylla_colors, pt.cex=2, cex=.8, bty="n", ncol=4)


plot(tensor.f.mb[[46]][["igraph"]], layout=tensor.f.mb[[46]][["coord"]],rescale=T, 
     vertex.size=tensor.f.mb[[46]][["vsize"]], 
     vertex.color=family_colors[as.factor(V(tensor.f.mb[[46]][["igraph"]])$family)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family,
     vertex.shape=shape[as.factor(V(tensor.f.mb[[46]][["igraph"]])$phylla)],
     #edge.width=1.3,
     main=names(tensor.f.mb[46]))
legend(x=-1, y=-1.1, levels(as.factor(V(tensor.f.mb[[46]][["igraph"]])$family)), pch=21,
       col="#777777", pt.bg=family_colors, pt.cex=2, cex=.8, bty="n", ncol=4)



plot(tensor.f.mb[[46]][["igraph"]], layout=tensor.f.mb[[46]][["coord"]],rescale=T, 
     vertex.size=tensor.f.mb[[46]][["vsize"]], 
     vertex.color=genus_colors[as.factor(V(tensor.f.mb[[46]][["igraph"]])$genus)],
     vertex.label=NA,#V(tensor.f.mb[[46]][["igraph"]])$family, 
     #vertex.shape=V(tensor.f.mb[[46]][["igraph"]])$family,
     edge.curved=.15,
     edge.width=2,
     main=names(tensor.f.mb[46]))
legend(x=-1, y=-1.1, levels(as.factor(V(tensor.f.mb[[46]][["igraph"]])$genus)), pch=21,
       col="#777777", pt.bg=genus_colors, pt.cex=2, cex=.8, bty="n", ncol=4)



l <- layout.circle(tensor.f.mb[[46]][["igraph"]])
plot(tensor.f.mb[[46]][["igraph"]], rescale=T,vertex.size=rowMeans(clr(tensor.m[[46]], 1))+6, 
     vertex.color=palette,vertex.label=V(tensor.f.mb[[46]][["igraph"]])$phylla, main=names(tensor.f.mb[46]))




legend(c(levels(as.factor(as.vector(unlist(unname(taxclass.mc["phylla",])))))),pch=21,ncol = 1)
plot(tensor.f.mb[[47]][["igraph"]], layout=tensor.f.mb[[47]][["coord"]], vertex.size=tensor.f.mb[[47]][["vsize"]], vertex.label=NA, main=names(tensor.f.mb[47]))
plot(tensor.f.mb[[48]][["igraph"]], layout=tensor.f.mb[[48]][["coord"]], vertex.size=tensor.f.mb[[48]][["vsize"]], vertex.label=NA, main=names(tensor.f.mb[48]))
plot(tensor.f.mb[[49]][["igraph"]], layout=tensor.f.mb[[49]][["coord"]], vertex.size=tensor.f.mb[[49]][["vsize"]], vertex.label=NA, main=names(tensor.f.mb[49]))


library(NetIndices)

sort(degree.distribution(tensor.f.mb[[46]][["igraph"]]))
shortest_paths(tensor.f.mb[[46]][["igraph"]],1:130)
test.g <- tensor.f.mb[[46]][["igraph"]]
test.adj <- get.adjacency(tensor.f.mb[[46]][["igraph"]],sparse=F)
# The function output consists of 10 network properties.
# $N            #number of nodes
#$Ltot        #number of links
#$LD        #link density (average # of links per node)
#$C            #the connectance of the graph
# This function measures connectance as L/(N*(N-1)) where L is links, and N is nodes
# Connectance can also be calculated as L/(N^2)
test.p <- GenInd(test.adj)

# The degree of a node refers to the number of links associated with a node.
# Degree can be measured as the links going in ("in degree"), out ("out degree"), or both.
# The degree() function takes a graph input and gives the degree of specified nodes.
# With the argument "v=V(graph)" you tell the function to give the degree of all nodes in the graph,
# while the "mode" argument specifies in, out, or both.
test.deg <- igraph::degree(test.g,v=V(test.g),mode="all")

# Degree distribution is the cumulative frequency of nodes with a given degree
# this, like degree() can be specified as "in", "out", or "all"
test.degdist <- degree.distribution(test.g,cumulative=T,mode="all")

# Using the power.law.fit() function I can fit a power law to the degree distribution
test.power <- power.law.fit(test.degdist)
# The output of the power.law.fit() function tells me what the exponent of the power law is ($alpha)
# and the log-likelihood of the parameters used to fit the power law distribution ($logLik)
# Also, it performs a Kolmogov-Smirnov test to test whether the given degree distribution could have
# been drawn from the fitted power law distribution.
# The function thus gives me the test statistic ($KS.stat) and p-vaule ($KS.p) for that test
# Then I can plot the degree distribution
par(mfrow=c(1,1))
par(mar = c(5, 5, 1, 1))
plot(test.degdist,log="xy",
     ylim=c(.01,10),
     bg="black",pch=21,
     xlab="Degree",
     ylab="Cumulative Frequency")
# And the expected power law distribution
lines(1:8,(1:8)^((-test.power$alpha)+1))
# Graphs typically have a Poisson distribution (if they are random),
# power law (preferential attachment), or truncated power law (many real networks) degree distribution


# Clustering coefficient is the proportion of
# a nodes neighbors that can be reached by other neighbors
# in igraph this property is apparently called "transitivity"
transitivity(test.g)
# gives the clustering coefficient of the whole network
transitivity(test.g,type="local")
# gives the clustering coefficient of each node
# Betweenness is the number of shortest paths between two nodes that go through each node of interest
graph.betweenness <- igraph::betweenness(test.g,v=V(test.g))
graph.edge.betweenness <- igraph::edge.betweenness(test.g,e=E(test.g))
# Closeness refers to how connected a node is to its neighbors
graph.closeness <- igraph::closeness(test.g,vids=V(test.g))
# Clustering coefficient, betweenness, and closeness
# all describe the small world properties of the network.
# A network with small world properties is one in which
# it takes a relatively short path to get from one node to the next
# (e.g., six degrees of separation)




head(tensor.f.mb[[1]][["igraph"]])


dd.mb <- degree.distribution(tensor.f.mbig[[1]]["igraph"])

plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,.35), type='b',col="red")#, 
#ylab="Frequency", xlab="Degree", main="Degree Distributions")
legend("topright", c("MB"), col=c( "red"), pch=1, lty=1)



se.mb.dfg1mc <- spiec.easi(df.g1.mc$data,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=50))
ig.mb.dfg1mc <- graph.adjacency(se.mb.dfg1mc$refit, mode='undirected')
vsize.dfg1mc <- rowMeans(clr(df.g1.mc$data, 1))+6
am.coord.dfg1mc <- layout.fruchterman.reingold(ig.mb.dfg1mc)





