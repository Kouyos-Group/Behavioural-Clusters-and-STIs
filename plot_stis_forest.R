# ------ OBJECTIVE: generate forest plots for nsCAI/nsP/syphilis/STI prediction models ---------------------------
# ------ DEPENDS DIRECTLY on analyse_stis_regression.R -----------------------------------------------------------

rm(list=ls())

# run analyse_stis_regression.R to define regression models 
source('analyse_stis_regression.R')


# choose behaviour variable and times of interest
behaviour <- 'nscai'          # "nsp" or "nscai"
studystart <- '20010701'    # beginning of period used for clustering
studyend <- '20170501'      # cut-off date
endend <- '20200630'        # end of study period


# prepare data for the desired set of clusters
df <- regressionPrep(behaviour, studystart, studyend, endend)

# -------------------------- define outcomes (cave: order matters!) --------------------------------------

outcomes <- list(paste0(behaviour, '_after'), 'std_after', 'syph_after', 'num_std_after',  'num_syph_after')
names <- c(name, 'Nurse/physician-reported STIs', 'Lab-confirmed syphilis',
           'Number of nurse/physician-reported STIs', 'Number of lab-confirmed syphilis episodes')



# PREPAREFORESTPLOTS: generates coefficient tables as basis for forest plots
# input: outcome (dependent variable, e.g. "nscai_after"), behavioural variable ("nscai"/"nsp"), dat (data),
#        print.summary (if TRUE, prints summary of each regression model)
# output: table with regression coefficients for given outcome

prepareForestplots <- function(outcome, behaviour, dat, print.summary = FALSE){
  
  runRegression(outcome, behaviour, dat, print.summary = FALSE)
  models <- list(m.2, m.3, m.4, m.5)
  if (print.summary == TRUE){print(lapply(models, function(x) summary(x)))}
  
  # get df with coefficients
  model_names <- NULL
  for (m in models){model_names <- c(model_names, as.character(m$formula[3]))}
  model_names <- c(1:4, model_names[2:length(model_names)])
  
  
  # (U)nivariable (R)egression output
  UR <- NULL # dummy line needed for appending
  for (m in models){UR <- rbind(UR, summary(m)$coefficients[2:nrow(summary(m)$coefficients),])}
  UR <- data.frame(UR[, c(1, 2, 4)]) 
  colnames(UR) <- c('est', 'se', 'pval')
  rownames(UR) <- model_names
  UR$lower <- UR$est - 1.96*UR$se
  UR$upper <- UR$est + 1.96*UR$se
  UR[, c('est', 'upper', 'lower')] <- exp(UR[, c('est', 'upper', 'lower')])
  UR$oddsCI <- paste0(round(UR$est, 2), " [", round(UR$lower, 2), ", ", round(UR$upper, 2), "]")
  UR$significance <- cut(UR$pval, c(0, 0.001, 0.01, 0.05, 1), right=FALSE, labels=c('***', '**', '*', ''))
  UR <- UR[, c('est', 'lower', 'upper', 'oddsCI', 'pval', 'significance')]
  for (r in c(1, 6, 8)){UR <- InsertRow(UR, rep(NA, ncol(UR)), r); UR[r, c(1:3)] <- 1}
  UR$model <- rep('univariable', nrow(UR))
  rownames(UR) <- NULL
  
  
  # (M)ultivariable (R)egression output
  MR <- summary(m.17)$coefficients[, c(1, 2, 4)]
  MR <- data.frame(MR[-1,])
  colnames(MR) <- c('est', 'se', 'pval')
  MR$lower <- MR$est - 1.96*MR$se
  MR$upper <- MR$est + 1.96*MR$se
  MR[, c('est', 'upper', 'lower')] <- exp(MR[, c('est', 'upper', 'lower')])
  MR$oddsCI <- paste0(round(MR$est, 2), " [", round(MR$lower, 2), ", ", round(MR$upper, 2), "]")
  MR$significance <- cut(MR$pval, c(0, 0.001, 0.01, 0.05, 1), right=FALSE, labels=c('***', '**', '*', ''))
  MR <- MR[, c('est', 'lower', 'upper', 'oddsCI', 'pval', 'significance')]
  for (r in c(1, 6, 8)){MR <- InsertRow(MR, rep(NA, ncol(MR)), r); MR[r, c(1:3)] <- 1}
  MR$model <- rep('multivariable', nrow(MR))
  rownames(MR) <- NULL
  
  
  # (R)egression (M)odel
  RM <- rbind(UR, MR)
  
  RM$names <- rep(c('Behavioural clusters', '', '', '', '', paste('Last', name), '', 'Prior syphilis','', 'Age'), 2)
  RM$subnames <- rep(c('0', '1', '2', '3', '4', 'Yes', 'No', 'Yes', 'No', ''), 2)
  RM$fullnames <- rep(c('BC 0', 'BC 1', 'BC 2', 'BC 3', 'BC 4', paste(name, 'No'), paste(name, 'Yes'), 'Syph No', 'Syph Yes', 'Age'), 2)
  
  RM_m <- RM[, c('fullnames', 'est', 'lower', 'upper', 'oddsCI', 'model')]
  RM_m$outcome <- rep(outcome, nrow(RM))
  RM_m$names <- rep(names[which(outcomes == outcome)], nrow(RM))
  
  # modify df for count outcomes
  if (substr(outcome, 1, 3) == 'num'){
    RM_m <- RM_m[RM_m$fullnames != 'Syph No',]
    RM_m$fullnames[RM_m$fullnames == 'Syph Yes'] <- 'Syph number'
  }
  
  return(RM_m)
  
}



# --------------- compute for all outcomes ----------

l1 <- lapply(outcomes, function (x)  prepareForestplots(x, behaviour = behaviour, dat = df))
regr <- do.call(rbind.data.frame, l1)

regr$fullnames <- factor(regr$fullnames, levels = rev(c('BC 0', 'BC 1', 'BC 2', 'BC 3', 'BC 4', paste(name, 'No'), paste(name, 'Yes'), 'Syph No', 'Syph Yes', 'Syph number', 'Age')))
regr$outcome <- factor(regr$outcome, levels=c("std_after","syph_after", paste0(behaviour, "_after"), "num_std_after", "num_syph_after"))
regr$names <- factor(regr$names, levels=c(name, "Nurse/physician-reported STIs", "Lab-confirmed syphilis", 
                                          "Number of nurse/physician-reported STIs", "Number of lab-confirmed syphilis episodes"))
levels(regr$names) <- c("Nurse/physician-reported STIs\n", "Lab-confirmed syphilis\n", paste(name, '\n'), 
                        "Number of \nnurse/physician-reported STIs\n", "Number of \nlab-confirmed syphilis episodes\n")

# split into binary and count outcomes
regr.bin <-  regr[substr(regr$outcome, 1, 3) != 'num',]
regr.cnt <-  regr[substr(regr$outcome, 1, 3) == 'num',]



# -------------------- plot --------------------------

getForest <- function(data){
  
  mPalette <- economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)] 
  t <- theme_minimal() + theme(legend.position = 'bottom',
                               legend.spacing.x = unit(0.2, 'line'),
                               legend.key.size = unit(0.6, 'line'),
                               #legend.text = element_text(size = 7),
                               legend.title = element_blank(),
                               axis.line.x = element_line(colour = 'black'),
                               axis.ticks.x = element_line(size = 0.25),
                               axis.title.x = element_text(hjust = 0.5),
                               aspect.ratio = 2.25,
                               panel.grid = element_blank(),
                               plot.title = element_text(hjust = 0.5),
                               #strip.text = element_text(face = 'bold'),
                               plot.caption = element_text(hjust = 0.5),
                               panel.spacing = unit(5, "lines"))
  
  forest <- ggplot(data) + geom_hline(yintercept = 1, colour = 'grey', linetype = 'solid') + 
    geom_pointrange(aes(x = fullnames, y = est, ymin = upper, ymax = lower, colour = model), 
                    shape = 15, position = position_dodge(width = 0.9), alpha = 0.9) +
    geom_text(data = data, mapping = aes(x = fullnames, y = Inf, label = oddsCI, position = model), hjust = 0, 
              size = 3, position = position_dodge(width = 0.9)) +
    
    # shaping
    scale_y_continuous(trans='log2') +
    coord_flip(clip = 'off') + facet_wrap(~ names, scales = "free_x") + t +
    scale_color_manual(values = mPalette) + guides(colour = guide_legend(reverse = TRUE)) +
    
    # labels
    labs(caption = "Numbers to the right of bars represent OR [95% CI].",
         x = '', y = 'Odds ratio') + 
    scale_x_discrete(breaks = c('BC 0', 'BC 1', 'BC 2', 'BC 3', 'BC 4', paste(name, 'No'), paste(name, 'Yes'), 'Syph No', 'Syph Yes', 'Syph number', 
                                'Age'),
                     labels = c('Behavioural clusters               0', '1', '2', '3', '4', 
                                paste(name, 'before cut-off             No'), 'Yes', 
                                'Syphilis before cutoff            No','Yes', 'Number of syphilis episodes    \nbefore cut-off                            ', 
                                'Age per ten years                     '))
  
  
  return(forest)
}


forest.bin <- getForest(regr.bin) 
forest.cnt <- getForest(regr.cnt) 

try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_forest_bin_', name, '.png'), 
    width = 30, height = 15, units = 'cm', res = 500)
forest.bin
dev.off()

try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_forest_cnt_', name, '.png'), 
    width = 25, height = 16, units = 'cm', res = 500)
forest.cnt
dev.off()



