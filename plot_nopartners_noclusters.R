# --------------- OBJECTIVE: generate plots with BICs for different numbers of clusters' association with number 
#--------------------------- of partners -----------------------------------------------------------------------
# --------------- DIRECTLY DEPENDS on analyse_nopartners_noclusters.R ------------------------------------------

rm(list=ls())

source('analyse_nopartners_noclusters.R')


# choose behaviour variable of interest
behaviour <- 'nsp'                          # "nsp" or "nscai"
studystart <- '20100101'                # beginning of period used for clustering
endend <- '20200630'                    # end of study period

# prepare data for desired cluster set
noClustersPrep.nopartners(behaviour, studystart, endend)

# nsCAI models
if (behaviour == "nscai"){
  ics.nscai.nonpred <- getICsForClusters.nopartners(read.tree('tree_nscai_excl_20100101_20160616_zphi'), 7, 'nscai', 'number_sexpartners', 
                                                    'Association of nsCAI clusters with number of partners', BIC, predictive = F)
  ics.nscai.pred <- getICsForClusters.nopartners(read.tree('tree_nscai_excl_20100101_20150531_zphi'), 7, 'nscai', 'number_sexpartners',
                                                 'nsCAI clusters to predict number of partners', BIC, predictive = TRUE)
  dfs <- list(ics.nscai.nonpred, ics.nscai.pred); name <- 'nsCAI'
  }


# nsP models 
if (behaviour == "nsp"){
  ics.nsp.nonpred <- getICsForClusters.nopartners(read.tree('tree_nsp_excl_20100101_20160616_zphi'), 7, 'nsp', 'number_sexpartners', 
                                                  'Association of nsP clusters with number of partners', BIC, predictive = F)
  ics.nsp.pred <- getICsForClusters.nopartners(read.tree('tree_nsp_excl_20100101_20150531_zphi'), 7, 'nsp', 'number_sexpartners', 
                                               'nsP clusters to predict number of partners', BIC, predictive = TRUE)
  dfs <- list(ics.nsp.nonpred, ics.nsp.pred); name <- 'nsP'
  }


df <- do.call(rbind.data.frame, dfs)
df_m <- melt(df, id.vars = c('outcome', 'no_clusters'))
df_m$outcome <- factor(df_m$outcome, levels=c(paste("Association of", name, "clusters with number of partners"),
                                              paste(name, "clusters to predict number of partners")))


df_plrt <- df_m[df_m$variable == 'p_LRT',]
df_m <- df_m[df_m$variable != 'p_LRT',]
df_plrt$value <- cut(df_plrt$value, c(0, 0.001, 0.01, 0.05, 0.1, 1), right=FALSE, labels=c('***', '**', '*', '.', ''))
df_plrt$value[is.na(df_plrt$value)] <- ''

lim <- 8
mPalette <- economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)]
t <- theme_minimal() + theme(legend.position = 'right',
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.title = element_blank(),
                             axis.ticks = element_line(size = 0.25),
                             aspect.ratio = 0.5,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5))

p <- ggplot(df_m[df_m$no_clusters <= lim, ], 
            aes(x = no_clusters, y = value, colour = variable, fill = variable)) +
  geom_line(size = 1) +  geom_point() +
  facet_wrap(~ outcome, scales = 'free_y') + 
  
  # labels and shaping
  scale_color_manual(values = mPalette, 
                     labels = c('Clusters + age', 
                                paste(name, '+ age'), 
                                paste('Clusters +', name, '+ age'))) +
  scale_fill_manual(values = mPalette, 
                    labels = c('Clusters + age', 
                               paste(name, ' + age'), 
                               paste('Clusters +', name, '+ age'))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.25, 0.05), add=c(0,0))) +
  geom_text(data = df_plrt[df_plrt$no_clusters <= lim,], 
            mapping = aes(x = no_clusters, y = - Inf,  label = value), size = 3, colour = 'black', vjust = -1) +
  guides(fill = FALSE) +
  labs(y = 'BIC', x = 'Number of clusters',
       caption = "\nAsterisks mark p value of LRT of full model without clusters vs. with clusters. *** <0.001; ** <0.01; * <0.05; . <0.1") + t +
  
  # adds "axis lines" to non-border facets
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

p


try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_nopartners_noclusters_', name,'.png'), width = 28, height = 7, units = 'cm', res = 500)
p
dev.off()
