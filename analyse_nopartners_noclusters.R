# -------------- OBJECTIVE: explore different numbers of clusters for number of partners analysis -----------
# -------------- DEPENDS ON: data_prep_zphi.R for data ------------------------------------------------------


source('0_functions.R')
source('datagen.R')

library(lmtest)
library(gmodels)
library(egg)
library(ggpubr)
library(dendextend)
library(ape)
library(lme4)



# ---------------- define noClustersPrep.nopartners: prepares data for desired cluster set -----------------------

noClustersPrep.nopartners <- function(behaviour, studystart, endend){
  
  writeLines("Preparing data for number of partners analysis with different numbers of clusters")
  
  # generate cluster sets and dataset if not done yet
  studyends <- c('20150531', '20160616')
  for (studyend in studyends){
    
    # create required cluster set if needed 
    if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, "zphi", sep = '_')) == FALSE){
      getBehaviouralClusters(behaviour, 'excl', zphi.only = TRUE, studystart, studyend, 2)}
    
    # make sure data has been prepared
    if(file.exists(paste('data_zphi', studystart, studyend, endend, sep = '_')) == FALSE){
      dataPrepZPHI(studystart, studyend, endend)}
  }
}



# ------------------------------- define GetICSForClusters.nopertners -------------------------------------------

# GETICSFORCLUSTERS.nopartners: calculates information criteria for regression models with different numbers of clusters
# input: hc.phylo = phylo object (tree) of desired cluster set
#        k = maximum number of clusters analysed
#        behaviour = behavioural variable used for clustering ("nscai"/"nsp")
#        outcome = dependent variable (e.g. "nscai_after")
#        title = title of resulting plot
#        FUN = information criterion function (defaults to BIC)
#        predictive = TRUE if cut-off is May 2015, before recording of sexual partners
# output: df with BIC and p_LRT values for different number of clusters

getICsForClusters.nopartners <- function(hc.phylo, k, behaviour, outcome, title, FUN = BIC, predictive){
  
  writeLines(paste("running regression with different numbers of clusters for", outcome))
  
  if (predictive == T){studyend <- '20150531'}
  if (predictive == F){studyend <- '20160616'}
  
  tail <- read.csv(paste('data_zphi', studystart, studyend, endend, sep = '_'))
  
  # workaround to get 'cluster 0 ppl' (MSM who fulill clustering inclusion criteria with no nsCAI/nsP)
  zero.clu <- read.csv(paste('cluster', behaviour,'excl', studystart, studyend, 'zphi', sep = '_')); zeros <- zero.clu[zero.clu$cluster == 0, ]
  
  ic <- NULL
  
  for (i in c(1:k)){
    
    cat("-")
    
    # make clusters
    clu <- cutree(hc.phylo, k = i) 
    cluster17 <- data.frame('id' = as.numeric(names(clu)), 'cluster' = as.numeric(clu))
    cluster17 <- rbind(cluster17, zeros)
    
    
    # merge with syph and determine model fit
    dat <- merge(cluster17, tail, by = 'id')
    
    
    m1 <- glmer(dat[, outcome] ~ age + as.factor(cluster) + (1|id), data = dat, family="poisson")#; summary(m1)
    m2 <- glmer(dat[, outcome] ~ age + dat[, behaviour] + (1|id), data = dat, family="poisson")#; summary(m2)
    m3 <- glmer(dat[, outcome] ~ age + dat[, behaviour] + as.factor(cluster) + (1|id), data = dat, family="poisson", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))#; summary(m3)

    
    b1 <- FUN(m1)
    b2 <- FUN(m2)
    b3 <- FUN(m3)
    plrt <- lrtest(m3, m2)$`Pr(>Chisq)`[2]
    
    ic <- rbind(ic, c(i + 1, b1, b2, b3, plrt))
    ic <- as.data.frame(ic); colnames(ic) <- c('no_clusters', 'BIC1', 'BIC2', 'BIC3', 'p_LRT')
    
  }
  
  
  
  m1 <- glmer(dat[, outcome] ~ age + (1|id), data = dat, family="poisson")#; summary(m1)
  m2 <- glmer(dat[, outcome] ~ age + dat[, behaviour] + (1|id), data = dat, family="poisson")#; summary(m2)
  m3 <- glmer(dat[, outcome] ~ age + dat[, behaviour] + (1|id), data = dat, family="poisson")#; summary(m3) 
  
  
  b1 <- FUN(m1)
  b2 <- FUN(m2)
  b3 <- FUN(m3)
  plrt <- lrtest(m3, m2)$`Pr(>Chisq)`[2]
  
  ic <- rbind(c(1, b1, b2, b3, plrt), ic)
  ic <- cbind(ic, 'outcome' = rep(title, nrow(ic)))
  
  ic
  
}



