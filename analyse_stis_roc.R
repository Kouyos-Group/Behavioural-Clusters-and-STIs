# --------------------- OBJECTIVE: generate ROC curves for predictivity analysis ------------------------------------


source('0_functions.R')
source('datagen.R')

library(pROC) 
library(randomForest) 
library(textclean)



ROCPrep <- function(behaviour, studystart, studyend, endend){

  writeLines("Preparing data for ROC analysis")
  
  # create required cluster set if needed and plot nsCAI/nsP trends
  if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) == FALSE){
    getBehaviouralClusters(behaviour, 'excl', zphi.only = FALSE, studystart, studyend, 4)}
  if (behaviour == 'nscai'){name <<- 'nsCAI'}
  if (behaviour == 'nsp'){name <<- 'nsP'}
  clu <- read.csv(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) 
  p <- plotTrends(behaviour, 'excl',  studystart, studyend, 2001, 2020, FALSE)
  writeLines(paste("Plot", name , "trends in clusters")); print(p)
  
  # make sure data has been prepared and merge it with cluster membreship data
  if(file.exists(paste('data_shcs', studystart, studyend, endend, sep = '_')) == FALSE){
    dataPrepSHCS(studystart, studyend, endend)}
  tail <- read.csv(paste('data_shcs', studystart, studyend, endend, sep = '_'))
  df <- merge(clu, tail, by = 'id')
  df$age <- as.numeric(cut(df$age, seq(20,90,10), right=FALSE, labels=c(1:7)))
  for (i in colnames(df)){if (is.element(i, c('std_after', 'num_std_after')) == TRUE){df[, i][is.na(df[, i])] <- 0}}
  for (i in colnames(df)){df <- df[!is.na(df[, i]),]}
  for (i in colnames(df)){  # make sure all categorical values are treated factors
    if (is.element(i, c('age', 'registration_year', 'registration_age', 'cd4_before', 'num_syph_after', 'num_syph_before', 
                        'num_std_after')) == FALSE){df[, i] <- as.factor(df[, i])}} 
  
  # create dichotomised variable for repeated STIs and syphilis 
  df$rep_syph_before <- 0
  df$rep_syph_before[as.numeric(df$num_syph_before) >= 2] <- 1
  df$rep_syph_after <- 0
  df$rep_syph_after[as.numeric(df$num_syph_after) >= 2] <- 1
  df$rep_std_after <- 0
  df$rep_std_after[as.numeric(df$num_std_after) >= 2] <- 1
  
  df

}

# ------------------------------------------ define getROCValues -----------------------------------------------

# GETROCVALUES: calculates receiver operater characteristic values for prediction models
# input: outcome = dependent variable (e.g. "nscai_after")
#        behaviour = behavioural variable used for clustering ("nsp"/"nscai)
#        dat = data
#        fam = family for regression
# output: df with ROC values and df with AUC values

getROCValues <- function(outcome, behaviour, dat, fam = 'binomial'){
  
  
  # univariable models
  m.1 <- glm(dat[, outcome] ~ NULL, family = fam, data = dat)
  m.2 <- glm(dat[, outcome] ~ cluster, family = fam, data = dat)
  m.3 <- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')], family = fam, data = dat)
  m.4 <- glm(dat[, outcome] ~ age, family = fam, data = dat)
  m.5 <- glm(dat[, outcome] ~ dat[,paste0(if (substr(outcome, 1, 3) != 'rep'){'num_'}, 'syph_before')], family = fam, data = dat)
  
  # multivariable models
  m.11 <- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')] + age, family = fam, data = dat)
  m.12 <- glm(dat[, outcome] ~ dat[,paste0(if (substr(outcome, 1, 3) == 'rep'){'rep_'}, 'syph_before')] + age, family = fam, data = dat)
  m.13 <- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')] + dat[,paste0(if (substr(outcome, 1, 3) == 'rep'){'rep_'}, 'syph_before')] + age, family = fam, data = dat)
  m.14 <- glm(dat[, outcome] ~ cluster + age, family = fam, data = dat)
  m.15 <- glm(dat[, outcome] ~ cluster + dat[, paste0(behaviour, '_before')] + age, family = fam, data = dat)
  m.16 <- glm(dat[, outcome] ~ cluster + dat[,paste0(if (substr(outcome, 1, 3) == 'rep'){'rep_'}, 'syph_before')] + age, family = fam, data = dat)
  m.17 <- glm(dat[, outcome] ~ cluster + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (substr(outcome, 1, 3) == 'rep'){'rep_'}, 'syph_before')] + age, family = fam, data = dat)
  
  
  models <- list(m.4, m.14, m.11, m.15, m.12, m.16, m.13, m.17)
  
  
  df.roc <- NULL
  df.auc <- NULL
  auc.text <- NULL
  for (m in models){
    out <- roc(dat[, outcome], m$fitted.values, plot=FALSE, ci = TRUE)
    roc <- data.frame(cbind("sens" = out$sensitivities, "spec" = out$specificities))
    roc <- cbind(rep(as.character(m$formula[3]), nrow(roc)), roc); colnames(roc)[1] <- 'model'
    roc <- cbind(rep(outcome, nrow(roc)), roc); colnames(roc)[1] <- 'outcome'
    df.roc <- rbind(df.roc, roc)
    
    auc <- data.frame(cbind("auc" = out$auc,  "lower" = out$ci[1], "upper" = out$ci[3]))
    auc <- cbind(rep(as.character(m$formula[3]), nrow(auc)), auc); colnames(auc)[1] <- 'model'
    auc <- cbind(rep(outcome, nrow(auc)), auc); colnames(auc)[1] <- 'outcome'
    df.auc <- rbind(df.auc, auc)
    
  }
  
  return(list(df.roc, df.auc))
  
}


# for arranging results and plotting refer to plot_stis_roc
