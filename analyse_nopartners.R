# -------------- OBJECTIVE: test whether clusters are associated with number of sexual partners ---------------------


source('0_functions.R')
source('datagen.R')


library(lme4)
library(lmtest)
library(stats)
library(utils)

regressionPrep.nopartners <- function(behaviour, studystart, studyend, endend, zphi.only){
  
  writeLines("-------------------------------------------------------------")
  writeLines("Preparing data for number of sexual partners regression analysis with desired cluster set")
  
  # ---------------------------------- DATA PREPARATION ---------------------------------------------------
  
  # create required cluster set if needed 
  if(file.exists(paste('cluster', behaviour, 'excl', studystart, studyend, if(zphi.only == T){'zphi'}, sep = '_')) == FALSE){
    getBehaviouralClusters(behaviour, 'excl', zphi.only = zphi.only, studystart, studyend, 2)}
  
  # read in cluster csv and plot nsCAI/nsP trends in clusters
  if (behaviour == 'nscai'){name <<- 'nsCAI'}
  if (behaviour == 'nsp'){name <<- 'nsP'}
  clu <- read.csv(paste('cluster', behaviour, 'excl', studystart, studyend, if(zphi.only == T){'zphi'}, sep = '_')) 
  p <- plotTrends(behaviour, 'excl',  studystart, studyend, 2006, 2019, zphi.only = zphi.only)
  writeLines(paste("Plot", name, "trends in clusters"))
  print(p)
  
  # make sure ZPHI data has been prepared and merge it with clusters 
  if(file.exists(paste('data_zphi', studystart, studyend, endend, sep = '_')) == FALSE){
    dataPrepZPHI(studystart, studyend, endend)}
  db <- read.csv(paste('data_zphi', studystart, studyend, endend, sep = '_'))
  db <- merge(clu, db, by = 'id')
  
  # exclude all with missing values 
  for (i in colnames(df)){df <- df[!is.na(df[, i]),]}
  db$age <- as.numeric(cut(db$age, seq(20,90,10), right=FALSE, labels=c(1:7)))
  for (i in colnames(db)){
    if (is.element(i, c('age', 'num_syph_before','number_sexpartners', 'mean_number_sexpartners')) == FALSE){
      db[, i] <- as.factor(db[, i])}} # make sure all non-continuous values are treated as factors 

  
  # print key figures
  writeLines("\n------------------------")
  writeLines("\nKEY FIGURES")
  writeLines(paste('Number of participants:', length(unique(db$id))))
  writeLines("Distribution of participants across clusters clusters [N (%)]:"); print(table(db$cluster, exclude = NULL)); print(round(prop.table(table(db$cluster, exclude = NULL))*100))
  writeLines("Summary of number of sexual partners"); print(summary(db$number_sexpartners))
  writeLines("Summary of number of sexual partners in cluster 0"); print(summary(db[db$cluster == 0,]$number_sexpartners))
  writeLines("Summary of number of sexual partners in cluster 1"); print(summary(db[db$cluster == 1,]$number_sexpartners))
  writeLines("Summary of number of sexual partners in cluster 2"); print(summary(db[db$cluster == 2,]$number_sexpartners))
  
  
  db
}

# ----------------------------------- define runRegression.nopartners -----------------------------------------------

# RUNREGRESSION: defines regression models to test whether clusters are associated with number of sexual partners
# input: behaviour = behavioural variable ("nscai"/"nsp")
#        print.summary = if TRUE, prints summary of each regression model (defaults to FALSE)
# output: list of regression models as glm objects

runRegression.nopartners <- function(behaviour, db, print.summary = FALSE){
  
  writeLines("\n------------------------")
  writeLines(paste("Running mixed effect Poisson and negative binomial regression for number of sexual partners with", 
                   name, "clusters inferred with data from", studystart, "until", studyend))
  
  # null model
  m0 <<- glmer(number_sexpartners ~ (1|id), data = db, family="poisson"); summary(m0)
  
  # poisson regression
  m1 <<- glmer(number_sexpartners ~ cluster + (1|id), data = db, family="poisson")
  m2 <<- glmer(number_sexpartners ~ db[, behaviour] + (1|id), data = db, family="poisson")
  m3 <<- glmer(number_sexpartners ~ age + (1|id), data = db, family="poisson")
  m4 <<- glmer(number_sexpartners ~ db[, behaviour] + age + (1|id), data = db, family="poisson")
  m5 <<- glmer(number_sexpartners ~ cluster + db[, behaviour] + age + (1|id), data = db, family="poisson", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  # negative binomial regression
  m6.1 <<- glmer(number_sexpartners ~ cluster + (1|id), data = db, family=MASS::negative.binomial(theta=10.0))
  m6.2 <<- glmer(number_sexpartners ~ db[, behaviour] + (1|id), data = db, family=MASS::negative.binomial(theta=10.0))
  m6.3 <<- glmer(number_sexpartners ~ db[, behaviour] + age + (1|id), data = db, family=MASS::negative.binomial(theta=10.0))
  m6.4 <<- glmer(number_sexpartners ~ cluster + db[, behaviour] + age + (1|id), data = db, family=MASS::negative.binomial(theta=10.0))
  
  models <<- list(m1, m2, m3, m4, m5, m6.1, m6.2, m6.3, m6.4)
  if (print.summary == TRUE){print(lapply(models, function(x) summary(x)))}
  
  writeLines("-------------------------------------------------------------")
  
}


# for arranging and plotting refer to plot_nopartners_forest.R

