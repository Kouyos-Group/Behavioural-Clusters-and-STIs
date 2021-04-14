# ---------------- OBJECTIVE: create and colour clusters dendrogram for methods figure --------------------


rm(list=ls())

source('0_functions.R')
source('datagen.R')


library(ggplot2)
library(reshape2)
library(ape)
library(phytools)
library(dendextend)
library(ggthemes)
library(lubridate)

behaviour <- 'nscai'
studystart <- '20010701'
studyend <- '20170501'
endend <- '20200630'


if (behaviour == 'nscai'){name <- 'nsCAI'}
if (behaviour == 'nsp'){name <- 'nsP'}

# create required cluster set if needed 
if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) == FALSE){
  getBehaviouralClusters(behaviour, 'excl', zphi.only = FALSE, studystart, studyend, 4)}
clusters.df <- read.csv(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) 
hc.phylo <- read.tree(paste('tree', behaviour, 'excl', studystart, studyend, sep = '_'))

# use to identify node labels: plot labelled dendrogram and draw an imaginary line through it vertically
# that intersects n edges (where n = number of clusters), then move right from the intersections and 
# write down number on next node - cave: slow!

#plot(hc.phylo, type =  "phylo", show.tip.label = FALSE, show.node.label = TRUE)
#nodelabels(cex = 0.5, bg = "white", horiz = TRUE)
#shownodelabs <- c((length(unique(clusters.df$id))+1):(length(unique(clusters.df$id))+15))

nodes <- c(5163, 4739, 3099, 2982)


# ------------------------------- colour clades (-> clusters) of tree -------------------------------------------




tree <- hc.phylo
for (n in nodes){
  tree <- paintSubTree(tree, node=n, state=as.character(which(n == nodes)+1), stem = TRUE)
}

mPalette <- c(economist_pal(fill = TRUE)(9)[c(2, 4, 5, 7, 8)]) #
names(mPalette) <- c(1:5)


try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_dendrogram_', name, '.png'), width = 10, height = 10, units = 'cm', res = 500)
plotSimmap(tree, mPalette, pts=F, lwd=1, fsize = 0.001, node.numbers=F, type = 'phylogram', direction = 'downwards') 
dev.off()



# plot nsCAI/nsP trends in clusters
p <- plotTrends(behaviour, 'excl',  studystart, studyend, 2001, 2020, FALSE); p
try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_cluster_trends_2001_2017_', name, '.png'), width = 10, height = 8, units = 'cm', res = 500)
p
dev.off()

