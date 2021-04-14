# ---- OBJECTIVE: generate bar plots with BIC (and p_LRT) values for models with different sets of predictors -------
# ---- DEPENDS DIRECTLY on analyse_stis_regression.R -----------------------------------------------------------------------------

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


# ------------------------------- define outcomes ---------------------------------------------

outcomes <- list(paste0(behaviour, '_after'), 'std_after', 'num_std_after', 'syph_after', 'num_syph_after')
names <- c(paste(name, 'after cut-off'), 'Nurse/physician-reported STIs after cut-off', 'Number of nurse/physician-reported STIs after cut-off', 
           'Lab-confirmed syphilis after cut-off', 'Number of lab-confirmed syphilis episodes after cut-off')



# COMPAREBICS: runs regression models and returns table with BIC and p_LRT values
# input: outcome (dependent variable, e.g. "nscai_after"), behavioural variable ("nscai"/"nsp"), dat (data),
#        print.summary (if TRUE, prints summary of each regression model)
# output: list with two dataframes - one with BIC values for each model and one with p_LRT for each 
#         with-clusters-vs-without-clusters comparison

compareBICs <- function(outcome, behaviour, dat, print.summary = FALSE){ 
  
  runRegression(outcome, behaviour, dat, print.summary = FALSE)
  models <- list(m.1, m.2, m.3, m.4, m.5, m.11, m.12, m.13, m.14, m.15, m.16, m.17)
  mv.models <- list(m.5, m.14, m.11, m.15, m.12, m.16, m.13, m.17)
  
  # likelihood ratio tests, BIC values
  bic <- unlist(lapply(mv.models, function(x) BIC(x)))
  p_lrt <- c(lrtest(m.5, m.14)$`Pr(>Chisq)`[2], lrtest(m.11, m.15)$`Pr(>Chisq)`[2], 
             lrtest(m.12, m.16)$`Pr(>Chisq)`[2], lrtest(m.13, m.17)$`Pr(>Chisq)`[2])
  
  if (outcome == 'nscai_after' | outcome == 'nsp_after'){
    namelis <<- unlist(lapply(mv.models, function(x) x$formula))  # somewhat messy but gives row names for tables
    clusters.inc <<- substr(namelis, 18, 24) == 'cluster'         # gives binary for whether model includes clusters
  }
  
  #print(outcome); print(bic)
  return(list(bic, p_lrt))
  
}



# ---------------------------------- generate table for BIC figure ------------------------------------ 

l3 <- lapply(outcomes, function (x)  compareBICs(x, behaviour = behaviour, dat = df))
bic <- do.call(rbind.data.frame, lapply(l3, `[[`, 1))
colnames(bic) <- substr(namelis, 18, nchar(namelis))
bic <- cbind(c(unlist(outcomes)), bic)


# make table with p_LRT for all model combinations (for annotating plot)
p_lrt <- do.call(rbind.data.frame, lapply(l3, `[[`, 2))
for(n in 1:nrow(p_lrt)){
  for(m in 1:ncol(p_lrt)){ifelse(as.numeric(p_lrt[n, m]) < 0.1, 
                                 p_lrt[n, m] <- as.character(cut(as.numeric(p_lrt[n,m]), 
                                                                 c(0, 0.001, 0.01, 0.05, 0.1), right=T, labels=c('p[LRT] <.001', 'p[LRT] <.01', 'p[LRT] <.05', 'p[LRT] <.1'))), #labels=c('<0.001', '<0.01', '<0.05', '<0.1'))),
                                 p_lrt[n, m] <- paste('p[LRT] ==', substr(as.character(round(as.numeric(p_lrt[n, m]), 3)), 2, 5)))
  }
}

colnames(p_lrt) <- c('Age', paste(name, 'before cut-off + age'), 'Syphilis before cut-off + age', paste(name, 'before cut-off + syphilis before cut-off + age'))
p_lrt <- cbind('outcome' = c(unlist(outcomes)), p_lrt)
p_lrt <- cbind(names, p_lrt)



# --------------------------------- organise results and plot -------------------------------------

# prepare p_lrt dataframe for annotation layer
p_lrt_m <- melt(p_lrt, id.vars = c('outcome', 'names'), variable.name = 'parameters')
p_lrt_m <- rbind(p_lrt_m, p_lrt_m)
p_lrt_m <- cbind(p_lrt_m, 'clusters.inc' = c(rep('without clusters', 20), rep('with clusters', 20)))


bic <- cbind(bic, names)
colnames(bic)[1] <- 'outcome'
bic_m <- melt(bic, id.vars = c('outcome', 'names'))
bic_m$clusters.inc <- "without clusters"
bic_m$clusters.inc[substr(as.character(bic_m$variable), 1, 7) == 'cluster'] <- "with clusters"
bic_m$parameters <- c(rep('Age', 2*nrow(bic)),
                      rep(paste(name, 'before cut-off + age'), 2*nrow(bic)), 
                      rep('Syphilis before cut-off + age', 2*nrow(bic)), 
                      rep(paste(name, 'before cut-off + syphilis before cut-off + age'), 2*nrow(bic)))


# make parameters and names factors to set order in which plotted
bic_m$parameters <- factor(bic_m$parameters, levels = rev(c('Age', 
                                                            paste(name, 'before cut-off + age'), 
                                                            'Syphilis before cut-off + age', 'Number of syphilis episodes before cut-off + age',
                                                            paste(name, 'before cut-off + syphilis before cut-off + age'), paste(name, 'before cut-off + number of syphilis episodes before cut-off + age'))))

bic_m$outcome <- factor(bic_m$outcome, levels=c("std_after","syph_after", paste0(behaviour, "_after"), "num_std_after", "num_syph_after"))
bic_m$names <- factor(bic_m$names, levels=c("Nurse/physician-reported STIs after cut-off", "Lab-confirmed syphilis after cut-off", paste(name, 'after cut-off'),
                                            "Number of nurse/physician-reported STIs after cut-off", "Number of lab-confirmed syphilis episodes after cut-off"))



# --------------------- plot figure -----------------------------

mPalette <- economist_pal(fill = TRUE)(9)
t <- theme_minimal() + theme(legend.position = c(0.83, 0.25),
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.title = element_blank(),
                             axis.ticks = element_line(size = 0.25),
                             aspect.ratio = 0.5,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5))

bic_plot <- ggplot(bic_m, aes(y = value, x = factor(parameters), fill = clusters.inc)) + 
  geom_bar(stat = "identity", position = 'dodge', width = 0.7, alpha = 0.95) +
  geom_text(data = p_lrt_m, mapping = aes(x = parameters, y = Inf, label = value), size = 2.5, hjust = 1, parse = TRUE) +
  
  # shaping
  coord_flip() + facet_wrap(~ names, scales = 'free_x') + t + guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = mPalette) + scale_x_discrete(expand = expand_scale(mult = c(.4, 0))) +
  scale_y_continuous(expand = expand_scale(mult = c(-0.6, 0.11), add=c(0,0))) +
  
  # labels
  labs(caption = paste('\nNumbers to right of bars represent p values of LRT between model with and without clusters. 
       For count outcomes (bottom row), syphilis before cut-off refers to number of syphilis episodes before cut-off.'),
       x = "Model parameters", y = 'BIC') + 
  
  # arrows and text for better/worse performance
  geom_segment(aes(x = 0.1, xend = 0.1, y = -Inf, yend = Inf), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  geom_segment(aes(x = 0.1, xend = 0.1, y = Inf, yend = -Inf), size = 0.25, arrow = arrow(length = unit(0.2, "line")), colour = 'grey') +
  annotate('text', x = 0.3, y = Inf, label = 'worse performance', size = 2.5, colour = 'grey', hjust = 1.1) +
  annotate('text', x = 0.3, y = -Inf, label = 'better performance', size = 2.5, colour = 'grey', hjust = -0.1) +
  
  # adds "axis lines" to non-border facets
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


#bic_plot

try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_bic_', name, '.png'), 
    width = 31, height = 13, units = 'cm', res = 500)
bic_plot
dev.off()
