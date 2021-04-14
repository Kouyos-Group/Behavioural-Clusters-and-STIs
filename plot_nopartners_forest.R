# ----- OBJECTIVE: run number of sexual partners regression and plot as forest plot --------------------
# ----- DEPENDS DIRECTLY on analyse_nopartners.R -------------------------------------------------------

rm(list=ls())

library(DataCombine) 
library(ggplot2) 
library(ggthemes) 
library(grDevices) 
library(lubridate) 


# run analyse_stis_regression.R to define regression models 
source('analyse_nopartners.R')


# choose behaviour variable and times of interest
behaviour <- 'nscai'              # "nsp" or "nscai"
studystart <- '20100101'          # beginning of period used for clustering
studyend <-  '20150531'           # cut-off date, '20160616' for overlapping analysis, '20150531' for predictive analysis
endend <- '20200630'              # end of study period
zphi.only <- TRUE                 # NULL if all ppl clustered, 'zphi' if only zphi ppl clustered


# prepare data for the desired set of clusters
db <- regressionPrep.nopartners(behaviour, studystart, studyend, endend, zphi.only)

# run regression models
runRegression.nopartners(behaviour, db)


# --------- results in a nutshell -----------

# LRT vs. null model
print(lrtest(m0, m1))              # LRT between model with clusters and null model

# LRT for poisson:
print(lrtest(m4, m5))              # LRT between model with clusters and model correcting for nsCAI/nsP and age

# LRTs for negative binomial:
print(lrtest(m6.3, m6.4))          # same with negative binomial model

# BICs: 
print(BIC(m4))                     # without clusters
print(BIC(m5))                     # with clusters





# ---------------- arrange output --------------------------

models <- list(m1, m2, m3)

# get df with coefficients
model_names <- NULL
for (m in models){model_names <- c(model_names, substr(as.character(summary(m)[15]), 38, nchar(as.character(summary(m)[15]))-43))}
model_names <- c(1:2, model_names[2:length(model_names)])


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
for (r in c(1, 4)){UR <- InsertRow(UR, rep(NA, ncol(UR)), r); UR[r, c(1:3)] <- 1}
UR$model <- rep('univariable', nrow(UR))
rownames(UR) <- NULL


# (M)ultivariable (R)egression output
MR <- summary(m5)$coefficients[, c(1, 2, 4)]
MR <- data.frame(MR[-1,])
colnames(MR) <- c('est', 'se', 'pval')
MR$lower <- MR$est - 1.96*MR$se
MR$upper <- MR$est + 1.96*MR$se
MR[, c('est', 'upper', 'lower')] <- exp(MR[, c('est', 'upper', 'lower')])
MR$oddsCI <- paste0(round(MR$est, 2), " [", round(MR$lower, 2), ", ", round(MR$upper, 2), "]")
MR$significance <- cut(MR$pval, c(0, 0.001, 0.01, 0.05, 1), right=FALSE, labels=c('***', '**', '*', ''))
MR <- MR[, c('est', 'lower', 'upper', 'oddsCI', 'pval', 'significance')]
for (r in c(1, 4)){MR <- InsertRow(MR, rep(NA, ncol(MR)), r); MR[r, c(1:3)] <- 1}
MR$model <- rep('multivariable', nrow(MR))
rownames(MR) <- NULL


# (R)egression (M)odel
RM <- rbind(UR, MR)


RM$names <- rep(c('Behavioural clusters', '', '',  paste('Last', name), '', 'Age'), 2)
RM$subnames <- rep(c('0', '1', '2', 'No', 'Yes',  ''), 2)
RM$fullnames <- rep(c('BC 0', 'BC 1', 'BC 2',  paste(name, 'No'), paste(name, 'Yes'),  'Age'), 2)
RM$fullnames <- factor(RM$fullnames, levels = rev(c('BC 0', 'BC 1', 'BC 2', paste(name, 'No'), paste(name, 'Yes'), 'Age')))

RM_m <- RM[, c('fullnames', 'est', 'lower', 'upper', 'oddsCI', 'model')]


mPalette <- economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)] 

t <- theme_minimal() + theme(legend.position = 'bottom',
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.title = element_blank(),
                             axis.title = element_text(size = 9),
                             axis.line.x = element_line(colour = 'black'),
                             axis.ticks.x = element_line(size = 0.25),
                             axis.title.x = element_text(hjust = 0.5),
                             aspect.ratio = 2.25,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'),
                             plot.caption = element_text(hjust = 0.5),
                             panel.spacing = unit(5, "lines"))

forest <- ggplot(RM_m) + geom_hline(yintercept = 1, colour = 'grey', linetype = 'solid') + 
  geom_pointrange(aes(x = fullnames, y = est, ymin = upper, ymax = lower, colour = model), 
                  shape = 15, position = position_dodge(width = 0.9), alpha = 0.9) +
  geom_text(data = RM_m, mapping = aes(x = fullnames, y = Inf, label = oddsCI, position = model), hjust = 0, #hjust = 1, 
            size = 3, position = position_dodge(width = 0.9)) +
  
  # shaping
  t +
  scale_y_continuous(trans='log2') + 
  coord_flip(clip = 'off') + 
  scale_color_manual(values = mPalette) + guides(colour = guide_legend(reverse = TRUE)) +
  
  # labels
  labs(caption = "Numbers to the right of bars represent IRR [95% CI].",
       title = paste('Factors associated with number of sexual partners\n', if (studyend == "20150531"){'after cut-off'}), 
       x = '', y = 'Incidence Rate Ratio \n(Exponent of coefficient in Poisson regression)') + 
  scale_x_discrete(breaks = c('BC 0', 'BC 1', 'BC 2', paste(name, 'No'), paste(name, 'Yes'), 'Age'),
                   labels = c('Behavioural clusters               0', '1', '2',
                              paste(name, 'before cut-off             No'), 'Yes',
                              'Age per ten years                     '),
                   expand = expand_scale(mult = c(.2, 0))) +
  
  # arrows and text for more/fewer partners
  geom_segment(aes(x = 0.1, xend = 0.1, y = 0.5, yend = Inf), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  geom_segment(aes(x = 0.1, xend = 0.1, y = 4, yend = 0), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  annotate('text', x = 0.3, y = Inf, label = 'more partners', size = 2.5, colour = 'grey', hjust = 1) +
  annotate('text', x = 0.3, y = 0, label = 'fewer partners', size = 2.5, colour = 'grey', hjust = 0) 



try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), "_nopartners_forest_", name, '_', studyend,".png"), width = 20, height = 14, units = 'cm', res = 500) 
forest
dev.off()


