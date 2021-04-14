# --------------- OBJECTIVE: explore predictive model with different numbers of clusters --------------------------- 


source('0_functions.R')
source('datagen.R')

library(lmtest)
library(gmodels)
library(egg)
library(ggpubr)
library(dendextend)
library(ape)


# ---------------- define noClustersPrep.nopartners: prepares data for desired cluster set -----------------------

noClustersPrep <- function(behaviour, studystart, studyend, endend){
  
  writeLines("Preparing data for nsCAI/nsP/STI/syphilis analysis with different numbers of clusters")
  
  # create required cluster set - or in this case mainly phylo object - if needed 
  if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) == FALSE){
    getBehaviouralClusters(behaviour, 'excl', zphi.only = FALSE, studystart, studyend, 4)}
  if (behaviour == 'nscai'){name <- 'nsCAI'}
  if (behaviour == 'nsp'){name <- 'nsP'}
  
  
  # make sure data has been prepared
  if(file.exists(paste('data_shcs', studystart, studyend, endend, sep = '_')) == FALSE){
    dataPrepSHCS(studystart, studyend, endend)}
  
  #  read in and prepare data and merge with cluster membership
  df <- read.csv(paste('data_shcs',studystart, studyend, endend, sep = '_'))
  df$age <- as.numeric(cut(df$age, seq(20,90,10), right=FALSE, labels=c(1:7)))
  for (i in colnames(df)){if (is.element(i, c('std_after', 'num_std_after')) == TRUE){df[, i][is.na(df[, i])] <- 0}}
  for (i in colnames(df)){df <- df[!is.na(df[, i]),]}
  
  # somewhat dirty workaround to get 'cluster 0 ppl' (MSM who fulill clustering inclusion criteria with no nsCAI/nsP)
  clu <- read.csv(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) 
  zeros <<- clu[clu$cluster == 0, ]
  
  df
  
}





# ------------------------------- define getICsForClusters -------------------------------------------

# GETICSFORCLUSTERS: calculates information criteria for regression models with different numbers of clusters
# input: hc.phylo = phylo object (tree) of desired cluster set
#        k = maximum number of clusters analysed
#        behaviour = behavioural variable used for clustering ("nscai"/"nsp")
#        fam = family for regression ("Poisson"/"binomial")
#        outcome = dependent variable (e.g. "nscai_after")
#        title = title of resulting plot
#        FUN = information criterion function (defaults to BIC)
# output: df with BIC and p_LRT values for different number of clusters


getICsForClusters <- function(hc.phylo, k, behaviour, outcome, fam, title, FUN = BIC){
  
  writeLines(paste("\nrunning", fam, "regression with different numbers of clusters for", outcome))
  ic <- NULL
  
  for (i in c(1:k)){
    
    cat("-")
    
    # make clusters
    clu <- cutree(hc.phylo, k = i) 
    cluster17 <- data.frame('id' = as.numeric(names(clu)), 'cluster' = as.numeric(clu))
    cluster17 <- rbind(cluster17, zeros)
    
    # merge with syph and determine model fit
    dat <- merge(cluster17, df, by = 'id')
    #if (syph.corr == TRUE){dat <- dat[!is.na(dat$ever_syph),]}
    
    #m0 <- glm(dat[, outcome] ~ age, family = fam, data = dat, na.action = 'na.omit')
    m1 <- glm(dat[, outcome] ~ age + as.factor(cluster), family = fam, data = dat, na.action = 'na.omit')
    m2 <- glm(dat[, outcome] ~ age + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')], family = fam, data = dat, na.action = 'na.omit')
    m3 <- glm(dat[, outcome] ~ age + as.factor(cluster) + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')], family = fam, data = dat, na.action = 'na.omit')
    
    b1 <- FUN(m1)
    b2 <- FUN(m2)
    b3 <- FUN(m3)
    plrt <- lrtest(m3, m2)$`Pr(>Chisq)`[2]
    
    ic <- rbind(ic, c(i + 1, b1, b2, b3, plrt))
    ic <- as.data.frame(ic); colnames(ic) <- c('no_clusters', 'BIC1', 'BIC2', 'BIC3', 'p_LRT')
  }
  
  
  # make NULL model line
  m1 <- glm(dat[, outcome] ~ age, family = fam, data = dat, na.action = 'na.omit')
  m2 <- glm(dat[, outcome] ~ age + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')], family = fam, data = dat, na.action = 'na.omit')
  m3 <- glm(dat[, outcome] ~ age + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')], family = fam, data = dat, na.action = 'na.omit')
  
  b1 <- FUN(m1)
  b2 <- FUN(m2)
  b3 <- FUN(m3)
  plrt <- lrtest(m3, m2)$`Pr(>Chisq)`[2] 
  
  ic <- rbind(c(1, b1, b2, b3, plrt), ic)
  
  ic <- cbind(ic, 'outcome' = rep(title, nrow(ic)))
  
  
  ic
}


# for plotting and arranging refer to plot_stis_noclusters.R





