# ------ OBJECTIVE: plot ROC curves and AUC bar charts for nsCAI/nsP/syphilis/STI prediction ------------------
# ------ DEPENDS DIRECTLY on analyse_stis_roc.R ---------------------------------------------------------------

rm(list=ls())

source('analyse_stis_roc.R')

# choose behaviour variable and times of interest
behaviour <- 'nsp'              # "nsp" or "nscai"
studystart <- '20010701'    # beginning of period used for clustering
studyend <- '20170501'      # cut-off date
endend <- '20200630'        # end of study period

# generate data for desired cluster set
df <- ROCPrep(behaviour, studystart, studyend, endend)


# -------------------------- select & name outcomes and run function -------------------------------

# CAVE: order matters!
outcomes <- list('std_after', 'syph_after',  paste0(behaviour, '_after'), 'rep_std_after', 'rep_syph_after')
names <- c('Nurse/physician-reported STIs', 'Lab-confirmed syphilis',  paste(name, ''),
           'More than one \nnurse/physician-reported STI', 'More than one \nlab-confirmed syphilis episode')

l1 <- lapply(outcomes, function(x) getROCValues(x, behaviour, df))


# ------------------------------ arrange output for plotting ----------------------------------------


# change model to parameters and clusters.inc column
roc.vals <- do.call(rbind.data.frame, lapply(l1, dplyr::nth, 1)) # last argument selects 1st element of each item in l1
roc.vals$clusters.inc <- 'without clusters'
roc.vals$clusters.inc[substr(roc.vals$model, 1, 7) == 'cluster'] <- 'with clusters'

roc.vals$parameters <- as.factor(sub('cluster \\+ ', '', roc.vals$model))
roc.vals$parameters <- factor(roc.vals$parameters, levels(roc.vals$parameters)[c(1, 2, 4, 3)])
levels(roc.vals$parameters) <- c('Age', 
                                 paste(name, 'before cut-off + age'), 
                                 'Syphilis before cut-off + age',  
                                 paste(name, 'before cut-off + syphilis before cut-off + age'))

roc.vals$outcome <- factor(roc.vals$outcome, levels = unlist(outcomes))
roc.vals$names <- roc.vals$outcome
levels(roc.vals$names) <- names


# ---------------------------------------- plot ROC curves ---------------------------------------------------


mPalette <- economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)] 

t <- theme_minimal() + theme(legend.position = c(0.85, 0.25),
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.text = element_text(size = 7),
                             legend.title = element_blank(),
                             axis.ticks = element_line(size = 0.25),
                             aspect.ratio = 1,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5),
                             panel.border = element_rect(fill = NA))

rocplot <- ggplot(data = roc.vals, aes(x = 1 - spec, y = sens)) + 
  geom_line(aes(linetype = clusters.inc, colour = parameters), lwd = 0.7) + 
  geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
  
  # shaping
  facet_wrap(~ names) + t + guides(colour = guide_legend(order = 0), linetype = guide_legend(reverse = TRUE, order = 1)) +
  scale_colour_manual(values = mPalette) +
  
  # labels
  labs(caption = "\nFor repeated STI/syphilis outcomes (bottom row), syphilis before cut-off refers to more than one syphilis episode before cut-off.",
       x = "False positive rate", y = 'True positive rate') #+


rocplot



try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_roc_', name, '.png'), 
    width = 0.8*27, height = 0.8*20, units = 'cm', res = 500)
rocplot
dev.off()





# -------------------------- configure and plot AUC values -------------------------------


auc.vals <- do.call(rbind.data.frame, lapply(l1, dplyr::nth, 2))
auc.vals$clusters.inc <- 'without clusters'
auc.vals$clusters.inc[substr(auc.vals$model, 1, 7) == 'cluster'] <- 'with clusters'

auc.vals$parameters <- as.factor(sub('cluster \\+ ', '', auc.vals$model))
auc.vals$parameters <- factor(auc.vals$parameters, levels(auc.vals$parameters)[c(1, 2, 4, 3)]) # changes and reverses order
levels(auc.vals$parameters) <- c('Age',
                                 paste(name, 'before cut-off + age'),
                                 'Syphilis before cut-off + age',
                                 paste(name, 'before cut-off + syphilis before cut-off + age'))

auc.vals$parameters <- factor(auc.vals$parameters, levels = rev(c('Age', 
                                                                  paste(name, 'before cut-off + age'), 
                                                                  'Syphilis before cut-off + age',
                                                                  paste(name, 'before cut-off + syphilis before cut-off + age'))))



auc.vals$outcome <- factor(auc.vals$outcome, levels = unlist(outcomes))
auc.vals$names <- auc.vals$outcome
levels(auc.vals$names) <- names



mPalette <- rev(economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)])
t <- theme_minimal() + theme(legend.position = c(0.83, 0.25),
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.text = element_text(size = 6),
                             legend.title = element_blank(),
                             axis.ticks = element_line(size = 0.25),
                             aspect.ratio = 0.5,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5))

aucplot <- ggplot(auc.vals) + 
  geom_bar(aes(y = auc, x = factor(parameters), fill = parameters, alpha = clusters.inc), stat = "identity", position = 'dodge', width = 0.7) + 
  geom_errorbar(aes(x = factor(parameters), ymin = lower, ymax = upper, alpha = clusters.inc), position = position_dodge(0.7), width = .4) +
  
  # shaping
  coord_flip() + facet_wrap(~ names) + t + 
  guides(fill = guide_legend(reverse = TRUE, order = 0), alpha = guide_legend(reverse = TRUE, order = 1)) +
  scale_fill_manual(values = mPalette) + scale_x_discrete(expand = expand_scale(mult = c(.4, 0))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0), add=c(0,0))) +
  scale_alpha_manual(values =c(0.9, 0.5)) +  
  
  # labels
  labs(caption = "\nFor repeated STI/syphilis outcomes (bottom row), syphilis before cut-off refers to more than one syphilis episode before cut-off.",
       x = "Model parameters", y = paste('Area under the ROC curve')) + 
  geom_segment(aes(x = 0.1, xend = 0.1, y = -Inf, yend = Inf), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  geom_segment(aes(x = 0.1, xend = 0.1, y = Inf, yend = -Inf), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  annotate('text', x = 0.3, y = Inf, label = 'better performance', size = 2, colour = 'grey', hjust = 1.1) +
  annotate('text', x = 0.3, y = -Inf, label = 'worse performance', size = 2, colour = 'grey', hjust = -0.1) +
  
  # adds "axis lines" to non-border facets
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


aucplot

try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_auc_', name, '.png'), 
    width = 27, height = 11, units = 'cm', res = 500)
aucplot
dev.off()

