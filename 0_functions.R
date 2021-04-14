# ----------------------------- collection of possibly handy functions -------------------------------------


# -------------------------------------- CLUSTERING FUNCTION -----------------------------------------------

# GETBEHAVIOURALCLUSTERS: cluster participants with hclust by given variable into n number of clusters
# input: behaviour = behavioural variable of interest ("nscai"/"nsp"), 
#        group = whether or not to exclude participants who never reported nsCAI/nsP during follow-up ("excl"/"incl"), 
#        zphi.only = whether or not to include only ZPHI participants in clustering algorithm (TRUE/FALSE; defaults to FALSE), 
#        startdate = startdate of period to be used for clustering as number string e.g. '20010107', 
#        enddate = enddate of period to be used for clustering as number string, 
#        n = number of clusters desired
# output: CSV with cluster classification and cluster dendrogram as phylo object
# note: CSV has clusters labelled 1-n, or 0-n if group = 'excl'


getBehaviouralClusters <- function(behaviour, group, zphi.only = F, startdate, enddate, n){
  
  writeLines(paste("Generating", behaviour, "clusters with data from", startdate, "to", enddate))
  
  source('datagen.R')
  
  library(readstata13) 
  library(tidyverse)
  library(plyr)
  library(graphics)
  library(gmodels)
  library(lubridate)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ape)
  library(phytools)
  library(dendextend)
  library(RColorBrewer)
  
  # functions:
  # creates trajectory for each patient, flipping in the middle of the cuts if the info changes
  
  pw.eval.ext <- function(cuts, data, tol = 10e-10, relaxation.period = 0.5){
    ## "data" is of class matrix
    ## "relaxation.period" is to extend the validity of records on the extremes
    
    #if only one registry, duplicate it
    if( dim(data)[1] == 1 ) data = rbind(data, data)
    
    #just that it looks like a dataframe
    data =  data.frame(t = data[,1], value = data[,2])
    
    #in case the registries aren't in chronological order:
    data = data[order(data$t),]
    
    print("ordered original data"); print(data)
    
    ##repeat the last regristry, "not to lose the information"
    data = rbind(data, c( data[nrow(data),1] + relaxation.period, data[nrow(data),2] ) )
    
    ##t to mid-points between registries + an auxiliar point in the past, at the midpoint between first registry and a hypotethical one, one relaxation.period before (records are retrospective by construction): (t's)
    data$t = c(max(data$t[1]-relaxation.period/2,0) ,(data$t[1:(length(data$t)-1)] + data$t[2:length(data$t)])/2)
    
    print("transformed data"); print(data)
    
    out = unlist(lapply(cuts, function(cut) {
      cut = cut + tol
      #locate the cut with respect to t_i and t_i+1 (intervals are of the form: t'0-t'1 -> v1 ... t'_i-t'_(i+1) -> v_(i+1)... t'_(n-1)-t'_(n) -> v_(n), we are looking for the relevant index to pick taht v
      ii.pick.value = which((data$t[c(1:(nrow(data)-1))] - rep(cut, nrow(data)-1))*(data$t[c(2:nrow(data))] - rep(cut, nrow(data)-1) ) <= 0.0)[1]
      if (cut < max(data$t)) data$value[ii.pick.value]
      else NA   ##beyond the data + extention, can't evaluate
    }
    ))
    out
  }
  
  
  # puts the by objects together as a data frame 
  by2dataframe.f =  function(by.ob) {for(i in 1:length(by.ob)){
    if(i == 1) {by.ob.df = by.ob[[1]]}else{by.ob.df = rbind(by.ob.df, by.ob[[i]])}}; by.ob.df}
  
  
  # ------------------------------ data pre-processing -------------------------------------------------
  
  if (file.exists("datagen") == F){
    generateData()
  }
  fup <- read.csv('datagen') # table from datagen with nscai/nsp from 2000 onwards
  fup$fupdate <- as.Date(fup$fupdate)
  fup$date <- as.numeric(gsub('-', '', fup$fupdate))
  fup <- fup[fup$date >= startdate & fup$date <= enddate,]
  fup <- fup[, c('id', 'fupdate', 'date', 'nsp', 'nscai')] 
  
  writeLines(paste('Number of participants in read-in for clustering:', length(unique(fup$id))))
  
  
  # ------------------------------- participant selection ------------------------------------------------
  
  
  # exclude ppl with <2 years of follow-up using follow-up matrix
  cuts = seq(as.numeric(min(fup$fupdate)), as.numeric(max(fup$fupdate)), by = 365.25/12) # x axis of fup matrix
  trajectories.follow.up.by = by(fup, fup$id, 
                                 function(m) pw.eval.ext(cuts, data = cbind(m$fupdate, 1)))
  
  trajectories.follow.up.df = data.frame(by2dataframe.f(trajectories.follow.up.by))
  rownames(trajectories.follow.up.df) <- unique(fup$id)
  
  # exclude ppl with fewer than 2 years of fup with valid nsp answer
  fup <- fup[is.element(fup$id, 
                        as.numeric(rownames(trajectories.follow.up.df
                                            [rowSums(!is.na(trajectories.follow.up.df)) >= 24, ]))),]
  writeLines(paste('Number of participants after excluding those with <2 years of fup:', length(unique(fup$id))))
  
  # exclude ppl with fewer than 2 valid nscai answers
  times_fup <- by(fup$nscai, fup$id, function(x) sum(!is.na(x), na.rm = T))
  times_fup <- data.frame('id' = as.numeric(names(times_fup)), 'times_fup' = as.numeric(times_fup))
  fup <- fup[is.element(fup$id, times_fup$id[times_fup$times_fup >= 2]),]
  writeLines(paste('Number of participants after excluding those with <2 valid condom use reports:', length(unique(fup$id))))
  
  
  # select only ZPHI participants if applicable 
  if (zphi.only == T){
    zphi <- read.csv('datagen_zphi')
    fup <- fup[is.element(fup$id, zphi$id), ]; writeLines('ZPHI participants selected')
    zphi <- 'zphi' # needed for naming of csvs etc.
  }
  if (zphi.only == F){zphi <- NULL}
  
  fup <- fup[!is.na(fup$id),]
  fup_original <- fup; length(unique(fup_original$id))
  
  
  # exclude those with no nscai reports if applicable 
  ev <- by(fup[, behaviour], fup$id, function(x) max(x, na.rm = T))
  ev <- data.frame(id = as.numeric(names(ev)), ever = as.numeric(ev))
  colnames(ev) <- c('id', paste('ever'))
  if (group == 'excl'){fup <- fup[is.element(fup$id, ev$id[ev$ever == 1]), ]}
  
  
  # -------------- create behavioural matrix and calculate pairwise distances between trajectories ------------- 
  
  
  writeLines(paste('Number of participants included in hclust:', length(unique(fup$id))))
  
  trajectories.per.id.by =  by(fup,
                               fup$id,
                               function(pat){
                                 data =  data.frame(as.numeric(pat$fupdate), pat[, as.character(behaviour)])
                                 out = pw.eval.ext(cuts, data = data)
                               })
  
  trajectories.per.id.matrix = by2dataframe.f(trajectories.per.id.by)
  rownames(trajectories.per.id.matrix) <- unique(fup$id)
  
  # calculate distance
  d0 =  dist(trajectories.per.id.matrix, method="binary")
  d0[is.na(d0)] = 1
  d0.m <- as.matrix(d0); nrow(d0.m)
  
  
  # ----------------------------------- run clustering algorithm -------------------------------------------
  
  
  hc = hclust(d0, method = "ward.D")
  
  # write hc to directory (as phylogenetic tree)
  hc.phylo <- as.phylo(hc)
  cluster.set.name <- paste(behaviour, group, startdate, enddate, sep = '_')
  if (zphi.only == TRUE){cluster.set.name <- paste(behaviour, group, startdate, enddate, 'zphi', sep = '_')}
  write.tree(hc.phylo, paste('tree', cluster.set.name, sep = '_'))
  writeLines(paste('Tree written: tree', cluster.set.name, sep = '_'))
  
  
  # ---------------------------- create dataframe with id and assigned cluster ----------------------------
  
  
  clu <- cutree(hc.phylo, k = n) 
  cluster_df <- data.frame('id' = as.numeric(names(clu)), 'cluster' = as.numeric(clu))
  
  if (group == 'excl'){
    fup_id <- data.frame(id=unique(fup_original[,'id']))  
    fup_id$cluster <- 0
    
    for (i in unique(fup$id)){
      if (i %in% unique(cluster_df$id)){
        fup_id$cluster[fup_id$id == i] <- cluster_df$cluster[cluster_df$id == i]}
    }
    cluster_df <- fup_id
  }
  
  
  # ------------------------- OUTPUT: csv file with cluster assignment  -------------------------------------
  
  writeLines(c(paste0('Total n = ', length(cluster_df$id), '. ', ncol(d0.m), ' participants clustered,', 
                     (length(cluster_df$id) - ncol(d0.m)), ' participants who never reported ', behaviour,' were added as cluster 0. '), 
               'The clusters contain the following proportions of the study population (rounded, in %): '))
  print(round(prop.table(table(cluster_df$cluster))*100, digits=1))
  
  # print table of nsCAI/nsP vs clusters as quality check
  # fupo <- merge(fup_original, cluster_df, by = 'id')
  # print(table(fupo$cluster, fupo[, behaviour]))
  
  write_csv(cluster_df, paste('cluster', cluster.set.name, sep = '_'))
  writeLines(c(paste('CSV written: cluster', cluster.set.name, sep = '_')))
  
  
  
  #cluster_df
}




# -------------------------------- nsCAI/nsP TREND PLOTTING FUNCTION ----------------------------------------

# PLOTTRENDS: calculates and plots nsCAI/nsP fraction for a given set of clusters
# input: behavioural variable of interest ("nscai" or "nsp"), whether or not to exclude participants who never 
#        reported nsCAI/nsP during follow-up ("excl"/"incl"), zphi.only (TRUE if only interested in ZPHI ppl), 
#        start- and enddate of period that was used for clustering as number string e.g. '20010107', 
#        start- and end year for plot (default 2001 to 2020)
# output: list of two ggplots, one for proportions and one for syph incidence

plotTrends <- function(behaviour, group, startdate, enddate, startyear = 2001, endyear = 2020, zphi.only = FALSE){
  
  
  library(readstata13) 
  library(tidyverse)
  library(plyr)
  library(graphics)
  library(gmodels)
  library(lubridate)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(grid)
  library(jtools)
  library(RColorBrewer)
  library(scales)
  
  
  variable <- paste(behaviour, group, sep ='_')
  if (zphi.only == TRUE){clu <- read.csv(paste('cluster', behaviour, group, startdate, enddate, 'zphi', sep = '_'))}
  if (zphi.only == FALSE){clu <- read.csv(paste('cluster', behaviour, group, startdate, enddate, sep = '_'))}
  if (behaviour == 'nscai'){name <- 'nsCAI'}
  if (behaviour == 'nsp'){name <- 'nsP'}
  
  minclu <- min(clu$cluster)
  maxclu <- max(clu$cluster)
  numclu <- maxclu - minclu + 1
  
  fup <- read.csv('datagen')
  fup$year <- format.Date(fup$fupdate, '%Y') 
  fup <- fup[fup$year >= 2000 & fup$year <= 2020,]
  names(fup)[names(fup) == behaviour] <- 'x'
  fup <- merge(fup, clu)
  length(unique(fup$id))
  
  
  # ------------------------- calculate nsCAI/nsP proportions over time -----------------------------------------
  
  
  # get NUMERATOR - number of [behaviour] reports per year
  by_year <- by(fup$x, fup$year, function(x) sum(x, na.rm=TRUE)) 
  by_year <- data.frame(year = as.numeric(names(by_year)), x_total = as.numeric(by_year))
  
  for (i in minclu:maxclu){
    clu_i <- getCount(fup, i, 'x')
    colnames(clu_i) <- c('year', paste('x', as.character(i), sep=''))
    by_year <- merge(by_year, clu_i, by='year', all = T)
  }
  
  
  # get DENOMINATOR - number of all reports per year
  total <- by(fup$timesfup, fup$year, sum) 
  total <- data.frame(year = as.numeric(names(total)), numtotal = as.numeric(total))
  by_year <- merge(by_year, total, by='year')
  
  for (i in minclu:maxclu){
    clu_i <- getCount(fup, i, 'timesfup')
    colnames(clu_i) <- c('year', paste('num', as.character(i), sep=''))
    by_year <- merge(by_year, clu_i, by='year')
  }
  
  # calculate [behaviour] proportion per year
  by_year$Total <- by_year$x_total/by_year$numtotal
  
  for (i in minclu:maxclu){
    by_year[as.character(i)] <- by_year[,which(c(minclu:maxclu) == i)+2]/by_year[,which(c(minclu:maxclu) == i)+3+numclu]
  }
  
  
  # ---------------------------------- calculate variance by year --------------------------------------------
  
  var_by_year <- data.frame('year' = unique(by_year$year))
  var_by_year$Total <- sqrt(by_year$Total*(1-by_year$Total)/by_year$numtotal)
  
  for (i in minclu:maxclu){
    var_by_year[as.character(i)] <- sqrt(by_year[,which(c(minclu:maxclu) == i)+4+2*numclu]*(1-by_year[,which(c(minclu:maxclu) == i)+4+2*numclu])/by_year[,which(c(minclu:maxclu) == i)+3+numclu])
  }
  
  
  by_year_m <- melt(by_year, 'year', c('Total', as.character(c(minclu:maxclu))), 
                    variable.name = 'cluster', value.name = 'x_proportion')
  var_by_year_m <- melt(var_by_year, 'year', c('Total', as.character(c(minclu:maxclu))), 
                        variable.name = 'cluster', value.name = 'variance')
  var_by_year_m$upper <- by_year_m[,3] + 1.96*var_by_year_m[,3]
  var_by_year_m$lower <- by_year_m[,3] - 1.96*var_by_year_m[,3]
  by_year_m <- cbind(by_year_m, upper = var_by_year_m[,4], lower = var_by_year_m[,5])
  
  
  
  # ------------------------------ plot proportion over time ----------------------------------------------
  
  mPalette <- c('lightgrey', economist_pal(fill = TRUE)(9)[c(2, 4, 5, 7, 8)])
  enddate.decimal = decimal_date(as.Date(paste(substr(enddate, 1, 4), substr(enddate, 5, 6), substr(enddate, 7, 8), sep = '-')))
  
  
  t <- theme(legend.position =  "right",
             legend.title =  element_blank(),
             legend.text =  element_text(size = 10),
             legend.key.size = unit(1,"line"),
             panel.background =  element_blank(),
             plot.title = element_text(size = 10, hjust = 0.5, face = 'bold'), #margin=margin(0,0,5,0),
             axis.text.x =  element_text(size = 10),
             axis.text.y =  element_text(size = 10),
             axis.title.x = element_text(size = 10, hjust = 0.5),#, margin=margin(10,0,0,0)),
             axis.title.y = element_text(size = 10),#, margin=margin(0,10,0,0)),
             axis.line.x =  element_line(colour = "black"),
             axis.line.y =  element_line(colour = "black"),
             aspect.ratio = 0.9
  )
  
  propplot <- ggplot(by_year_m[by_year_m$year >= startyear & by_year_m$year <= endyear,], 
                     aes(x =  year + 0.5, y = x_proportion, colour = as.factor(cluster), fill = as.factor(cluster))) +
    geom_vline(xintercept = seq(startyear + 1, endyear, by = 1), col  = "grey", lty = 2) +
    annotate(geom = "rect", xmin = enddate.decimal, xmax = endyear + 0.5, ymin = -2, ymax = 2, fill = "blue", alpha = 0.15) +
    geom_line(size=0.7) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA ) +
    t +
    scale_colour_manual(values = mPalette) +  scale_fill_manual(values = mPalette) + 
    ylab(paste("Fraction of men who report", name, sep = ' ')) + xlab("Calendar year")  + 
    ggtitle(paste(name, 'trends in behavioural clusters')) +
    coord_cartesian(xlim = c(startyear + 1, endyear + 0.5), ylim = c(0,1)) +
    geom_vline(xintercept = enddate.decimal, col  = "black", lty = 'longdash') 
  
  propplot  
}


# ---------------------------------- DATA PREPARATION FUNCTIONS -----------------------------------------

# GETPROPORTION: gets proportion of a variable over time
# input: dat = data, x = numerator, y = denominator
# output: dataframe with columns [year], numerator, denominator, proportion, upper and lower binomial error (95% CI) 
# (originally used in visualising sexual risk behaviour over time)

getProportion <- function(x, y = 'timesfup', dat = fup){
  print(x)
  dat <- dat[!is.na(dat[,x]),]
  numr <- by(dat[,x], dat$year, sum)
  numr <- data.frame(year = as.numeric(names(numr)), total = as.numeric(numr))
  colnames(numr) <- c('year', paste('numr', as.character(x), sep='_'))
  denr <- by(dat[,y], dat$year, sum)
  denr <- data.frame(year = as.numeric(names(denr)), total = as.numeric(denr))
  colnames(denr) <- c('year', paste('denr', as.character(x), sep='_'))
  prop <- merge(numr, denr, by='year')
  prop[as.character(x)] <- numr[,2]/denr[,2]
  prop[paste('upper', x, sep='_')] <- prop[as.character(x)] + 1.96 * sqrt(prop[as.character(x)]*(1-prop[as.character(x)])/prop[,3])
  prop[paste('lower', x, sep='_')] <- prop[as.character(x)] - 1.96 * sqrt(prop[as.character(x)]*(1-prop[as.character(x)])/prop[,3])
  prop[,c(2:6)]
}


# GETDENSITY: gets proportion of a variable per person (id)
# input: dat = data, x = numerator, y = denominator
# output: dataframe with columns id, numerator, denominator, proportion
# (originally used in regression for risk behaviour and syphilis)

getDensity <- function(x, y = 'timesfup', dat = fup){
  print(x)
  print(length(unique(dat$id)))
  numr <- by(dat[,x], dat$id, sum)
  numr <- data.frame(year = as.numeric(names(numr)), total = as.numeric(numr))
  colnames(numr) <- c('id', paste('numr', as.character(x), sep='_'))
  denr <- by(dat[,y], dat$id, sum)
  denr <- data.frame(year = as.numeric(names(denr)), total = as.numeric(denr))
  colnames(denr) <- c('id', paste('denr', as.character(x), sep='_'))
  dens <- merge(numr, denr, by='id')
  dens[paste0(as.character(x), '_den')] <- numr[,2]/denr[,2]
  dens
}


# GETEVER: gets a binary for whether a person has a hit for a variable (0 if never, 1 if ever)
# input: x = value to find ever of, dat = data
# output: dataframe with columns id and ever

getEver <- function(df, x, datename, studystart, studyend){
  #print(x)
  studyend <- as.Date(studyend)
  studystart <- as.Date(studystart)
  colnames(df)[colnames(df) == datename] <- 'date'
  df <- df[df$date <= studyend & df$date >= studystart,]
  ev <- by(df[,x], df$id, function(x) max(x, na.rm = TRUE))
  ev <- data.frame(id = as.numeric(names(ev)), ever = as.numeric(ev))
  ev$ever[ev$ever == -Inf] <- NA 
  colnames(ev) <- c('id', paste('ever', x, sep='_'))
  #print(length(ev$id))
  ev
}


# GETNUMEVENTS: gets a binary for whether a person has a hit for a variable (0 if never, 1 if ever)
# input: x = value to find ever of, dat = data
# output: dataframe with columns id and ever

getNumEvents <- function(df, x, datename, studystart, studyend){
  #print(x)
  studyend <- as.Date(studyend)
  studystart <- as.Date(studystart)
  colnames(df)[colnames(df) == datename] <- 'date'
  df <- df[df$date <= studyend & df$date >= studystart,]
  ev <- by(df[,x], df$id, function(x) sum(x, na.rm = TRUE))
  ev <- data.frame(id = as.numeric(names(ev)), ever = as.numeric(ev))
  ev$ever[ev$ever == -Inf] <- NA 
  colnames(ev) <- c('id', paste('num', x, sep='_'))
  #print(length(ev$id))
  ev
}


# GETCOUNT: gets count of a certain variable per year in given cluster (basically an extended by function)
# input: dat = data, i = cluster number, cvar = variable to be counted
# output: dataframe with columns year and count

getCount <- function(dat, i, cvar){
  dat_i <- subset(dat, cluster == i)
  by_year_i <- by(dat_i[, cvar], dat_i$year, function(x) sum(x, na.rm = TRUE))
  by_year_i <- data.frame(year = as.numeric(names(by_year_i)), num_events = as.numeric(by_year_i))
  by_year_i
}


# FINDLASTVALUE: gets last value per id for variable of interest up to date of interest
# input: dataframe, variable to get last of, name of date column, end of study period
# output: 2 * n(patients) df with id and last value (named e.g. 'cd4_2017)

findLastValue <- function(df, x, datename, studystart, studyend){
  studyend <- as.Date(studyend)
  studystart <- as.Date(studystart)
  colnames(df)[colnames(df) == datename] <- 'date'
  df <- df[df$date <= studyend & df$date >= studystart,]
  d <- df[!is.na(df[, x]),]
  x_last <- d[tapply(1:nrow(d), d$id, function(ii) ii[which.max(d$date[ii])]),]
  x_last$x_last <- x_last[, x]
  x_last <- x_last[, c('id', 'x_last')]
  colnames(x_last)[colnames(x_last) == 'x_last'] <- paste0(x, '_before') #paste(x, format.Date(studyend, '%Y'), sep = '_')
  x_last
}


# FINDFIRSTVALUE: gets first value per id for variable of interest up to date of interest
# input: dataframe, variable to get first of, name of date column, end of study period
# output: 2 * n(patients) df with id and first value (named e.g. 'cd4_2017')

findFirstValue <- function(df, x, datename, studystart, studyend){
  studyend <- as.Date(studyend)
  studystart <- as.Date(studystart)
  colnames(df)[colnames(df) == datename] <- 'date'
  df <- df[df$date <= studyend & df$date >= studystart,]
  d <- df[!is.na(df[, x]),]
  x_first <- d[tapply(1:nrow(d), d$id, function(ii) ii[which.min(d$date[ii])]),]
  x_first$x_first <- x_first[, x]
  x_first <- x_first[, c('id', 'x_first')]
  colnames(x_first)[colnames(x_first) == 'x_first'] <- paste0(x, '_after') #paste(x, format.Date(studystart, '%Y'), sep = '_')
  x_first
}


# PRINTNUMPROP: prints number and proportion of binary variable in desired format
# input: variable of interest (column name, as character), and datafram
# output: character with number of hits and proportion in brackets

printNumProp <- function(variable, data){
  t <- table(data[, variable], exclude = NULL)
  p <- round(prop.table(table(data[, variable], exclude = NULL))*100)
  paste0(t[[2]], " (", p[[2]], "%)")
}






