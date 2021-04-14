# -------------- OBJECTIVE: analyse whether clusters predict nsCAI, nsP, nurse/physician-reported STIs --------------
#------------------------------------- and laboratory-confirmed syphilis --------------------------------------------


source('0_functions.R')
source('datagen.R')


library(readstata13) 
library(lubridate)
library(lmtest)
library(stats)
library(utils)
library(DataCombine)
library(ggplot2)
library(ggthemes)
library(grDevices)
library(reshape2)


# ---------------- define regressionPrep (function that prepares data with desired cluster set) -------------------

regressionPrep <- function(behaviour, studystart, studyend, endend){
  
  writeLines("-------------------------------------------------------------")
  writeLines("Preparing data for nsCAI/nsP/STI/syphilis prediction analysis with desired cluster set")
  if (behaviour == 'nscai'){name <<- 'nsCAI'}
  if (behaviour == 'nsp'){name <<- 'nsP'}
  
  # make sure desired clusters have been inferred and plot nsCAI/nsP trends in clusters
  if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) == FALSE){
    getBehaviouralClusters(behaviour, 'excl', zphi.only = FALSE, studystart, studyend, 4)}
  clu <- read.csv(paste('cluster', behaviour, 'excl', studystart, studyend, sep = '_')) 
  p <- plotTrends(behaviour, 'excl',  studystart, studyend, 2001, 2020, FALSE)
  writeLines(paste("Plot", name , "trends in clusters")); print(p)
  
  # make sure SHCS data has been generated and add cluster membership, plus fine-prepare data
  if(file.exists(paste('data_shcs', studystart, studyend, endend, sep = '_')) == FALSE){
    dataPrepSHCS(studystart, studyend, endend)}
  tail <- read.csv(paste('data_shcs', studystart, studyend, endend, sep = '_'))
  df <- merge(clu, tail, by = 'id')
  df$age.original <- df$age; df$age <- as.numeric(cut(df$age, seq(20,90,10), right=FALSE, labels=c(1:7)))
  for (i in colnames(df)){if (is.element(i, c('std_after', 'num_std_after')) == TRUE){df[, i][is.na(df[, i])] <- 0}}
  for (i in colnames(df)){if (is.element(i, c('age', 'age.original', 'num_syph_after', 'num_syph_before', 'num_std_after')) == FALSE){df[, i] <- as.factor(df[, i])}} 
  for (i in colnames(df)){df <- df[!is.na(df[, i]),]}
  
  writeLines("\n------------------------")
  writeLines("\nPOPULATION CHARACTERISTICS:")
  writeLines(paste('Number of participants:', length(unique(df$id))))
  writeLines(paste('Median age:', median(df$age.original)))
  writeLines(paste("Reported nsP during follow-up:", printNumProp("ever_nsp", df)))
  writeLines(paste("Reported nsCAI during follow-up:", printNumProp("ever_nscai", df)))
  writeLines(paste("Laboratory-confirmed syphilis before cut-off:", printNumProp("syph_before", df)))
  writeLines(paste("More than one laboratory-confirmed syphilis episode before cut-off:", paste0(table(df$num_syph_before>1)[[2]], ' (', round(prop.table(table(df$num_syph_before>1))*100)[[2]], '%)')))
  writeLines("\nKEY FIGURES FOR ANALYSIS:")
  writeLines("Distribution of participants across clusters clusters [N (%)]:"); print(table(df$cluster, exclude = NULL)); print(round(prop.table(table(df$cluster, exclude = NULL))*100))
  writeLines(paste("Reported nsCAI after cut-off [N (%)]:", printNumProp("nscai_after", df)))
  writeLines(paste("Reported nsP after cut-off [N (%)]:",  printNumProp("nsp_after", df)))
  writeLines(paste("Nurse/physician-reported STI after cut-off [N (%)]:",  printNumProp("std_after", df)))
  writeLines(paste("Laboratory-confirmed syphilis after cut-off [N (%)]:",  printNumProp("syph_after", df)))
  writeLines(paste("More than one nurse/physician-reported STI after cut-off [N (%)]:",  paste0(table(df$num_std_after>1)[[2]], ' (', round(prop.table(table(df$num_std_after>1))*100)[[2]], '%)')))
  writeLines(paste("More than one laboratory-confirmed syphilis episode after cut-off [N (%)]:", paste0(table(df$num_syph_after>1)[[2]], ' (', round(prop.table(table(df$num_syph_after>1))*100)[[2]], '%)')))
  writeLines("-------------------------")
  
  df 
}



# ----------------------------------------- define runRegression -----------------------------------------------

# RUNREGRESSION: defines regression models to test whether clusters are predictive of nsCAI/nsP and STIs
# input: outcome = dependent variable (e.g. "nscai_after") 
#        behaviour = behavioural variable ("nscai"/"nsp")
#        dat = data
#        print.summary = if TRUE, prints summary of each regression model (defaults to FALSE)
# output: list of regression models as glm objects

runRegression <- function(outcome, behaviour, dat, print.summary = FALSE){
  
  # if outcome is number of STIs/syphilis episodes (i.e. "num_std_after"/"num_syph_after"), use Poisson regression
  if (substr(outcome, 1, 3) == 'num'){fam = 'poisson'}
  else {fam = 'binomial'}
  writeLines(paste("Running", fam, "regression for", outcome))
  
  # univariable models
  m.1 <<- glm(dat[, outcome] ~ NULL, family = fam, data = dat)
  m.2 <<- glm(dat[, outcome] ~ cluster, family = fam, data = dat)
  m.3 <<- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')], family = fam, data = dat)
  m.4 <<- glm(dat[, outcome] ~ dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')], family = fam, data = dat)
  m.5 <<- glm(dat[, outcome] ~ age, family = fam, data = dat)
  
  # multivariable models
  m.11 <<- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')] + age, family = fam, data = dat)
  m.12 <<- glm(dat[, outcome] ~ dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')] + age, family = fam, data = dat)
  m.13 <<- glm(dat[, outcome] ~ dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')] + age, family = fam, data = dat)
  m.14 <<- glm(dat[, outcome] ~ cluster + age, family = fam, data = dat)
  m.15 <<- glm(dat[, outcome] ~ cluster + dat[, paste0(behaviour, '_before')] + age, family = fam, data = dat)
  m.16 <<- glm(dat[, outcome] ~ cluster + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')] + age, family = fam, data = dat)
  m.17 <<- glm(dat[, outcome] ~ cluster + dat[, paste0(behaviour, '_before')] + dat[,paste0(if (fam == 'poisson'){'num_'}, 'syph_before')] + age, family = fam, data = dat)
  
  models <- list(m.1, m.2, m.3, m.4, m.5, m.11, m.12, m.13, m.14, m.15, m.16, m.17)
  if (print.summary == TRUE){print(lapply(models, function(x) summary(x)))}
  
  
}


# refer to plot_shcs_forestplot.R, table_plrt_ics.R and plot_bic_barplot.R for visualisations of regression results


