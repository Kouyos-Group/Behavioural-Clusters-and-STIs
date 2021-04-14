# ------ OBJECTIVE: summarise results from model comparison in table with p_LRT and AIC/BIC values ---------------
# ------ DEPENDS DIRECTLY on analyse_stis_regression.R ------------------------------------------------------------------------

rm(list=ls())

# run analyse_stis_regression.R to define regression models 
source('analyse_stis_regression.R')


# choose behaviour variable and times of interest
behaviour <- 'nscai'        # "nsp" or "nscai"
studystart <- '20010701'    # beginning of period used for clustering
studyend <- '20170501'      # cut-off date
endend <- '20200630'        # end of study period


# prepare data for the desired set of clusters
df <- regressionPrep(behaviour, studystart, studyend, endend)

# define outcomes 
outcomes <- list(paste0(behaviour, '_after'), 'std_after', 'num_std_after', 'syph_after', 'num_syph_after')
names <- c(paste(name, 'after cut-off'), 'Nurse/physician-reported STIs after cut-off', 'Number of nurse/physician-reported STIs after cut-off', 
           'Lab-confirmed syphilis after cut-off', 'Number of lab-confirmed syphilis episodes after cut-off')




# TESTCLUSTERCONTRIBUTION: runs regression models as specified and returns table with p_LRT/AIC/BIC values
# input: outcome (dependent variable, e.g. "nscai_after"), behavioural variable ("nscai"/"nsp"), dat (data),
#        print.summary (if TRUE, prints summary of each regression model)
# output: p_LRT and information criteria for given outcome (corresponds to row in Table 2)

testClusterContribution <- function(outcome, behaviour, dat, print.summary = FALSE){
  
  runRegression(outcome, behaviour, dat, print.summary = FALSE)
  models <- list(m.1, m.2, m.3, m.4, m.5, m.13, m.14, m.17)
  
  # likelihood ratio tests, information criteria
  ics <- c(lrtest(m.1, m.2)$`Pr(>Chisq)`[2], lrtest(m.13, m.17)$`Pr(>Chisq)`[2], round(c(BIC(m.13), BIC(m.17), AIC(m.13), AIC(m.17))))
  print(outcome); print(ics)
  return(ics)
}



# --------- generate information criteria table (Table 2) ---------
l2 <- lapply(outcomes, function (x)  testClusterContribution(x, behaviour = behaviour, dat = df))
ics <- do.call(rbind.data.frame, l2)
ics <- cbind(unlist(outcomes), ics)
colnames(ics) <- c('outcome', 'p_LRT_uv', 'p_LRT_mv', 
                   'BIC_without_clusters', 'BIC_with_clusters', 'AIC_without_clusters', 'AIC_with_clusters')

tbl <- ics
tbl$p_LRT_vs_null_model <- as.character(cut(tbl$p_LRT_uv, c(0, 0.001, 0.01, 0.05, 0.1), right=T, labels=c('<0.001', '<0.01', '<0.05', '<0.1')))
tbl$p_LRT_without_vs_with_clusters <- as.character(cut(tbl$p_LRT_mv, c(0, 0.001, 0.01, 0.05, 0.1), right=T, labels=c('<0.001', '<0.01', '<0.05', '<0.1')))
tbl$p_LRT_vs_null_model[is.na(tbl$p_LRT_vs_null_model)] <- round(tbl$p_LRT_uv[is.na(tbl$p_LRT_vs_null_model)], 3)
tbl$p_LRT_without_vs_with_clusters[is.na(tbl$p_LRT_without_vs_with_clusters)] <- round(tbl$p_LRT_mv[is.na(tbl$p_LRT_without_vs_with_clusters)], 3)
tbl$outcome <- c(name,'Any nurse/physician-reported STI', 'Number of nurse/physician-reported STIs',
                 'Laboratory-confirmed syphilis', 'Number of laboratory-confirmed syphilis episodes')
tbl <- tbl[, c('outcome', 'p_LRT_vs_null_model', 'p_LRT_without_vs_with_clusters', 
               "AIC_without_clusters", "AIC_with_clusters", "BIC_without_clusters",  "BIC_with_clusters")]
colnames(tbl) <- c('Outcome', 'p_LRT vs null model', 'p_LRT without vs with clusters', 
                   "AIC without clusters", "AIC with clusters", "BIC without clusters",  "BIC with clusters")

View(tbl)

