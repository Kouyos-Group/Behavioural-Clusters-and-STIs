# --------------- OBJECTIVE: generate plots with BICs for different numbers of clusters ------------------------
# --------------- DIRECTLY DEPENDS on noclusters_shcs.R! -------------------------------------------------------

rm(list=ls())

source('analyse_stis_noclusters.R')


# choose behaviour variable and times of interest
behaviour <- 'nsp'          # "nsp" or "nscai"
studystart <- '20010701'    # beginning of period used for clustering
studyend <- '20170501'      # cut-off date
endend <- '20200630'        # end of study period

# prepare data with desired cluster set
df <- noClustersPrep(behaviour, studystart, studyend, endend)


if (behaviour == 'nscai'){
  ics.nscai.nscai <- getICsForClusters(read.tree('tree_nscai_excl_20010701_20170501'), 7, 'nscai', 'nscai_after', 'binomial', 'nsCAI', BIC)
  ics.nscai.std <- getICsForClusters(read.tree('tree_nscai_excl_20010701_20170501'), 7,  'nscai','std_after', 'binomial', 'Nurse/physician-reported STIs', BIC)
  ics.nscai.num.std <- getICsForClusters(read.tree('tree_nscai_excl_20010701_20170501'), 7,  'nscai','num_std_after', 'poisson', 'Number of nurse/physician-reported STIs', BIC)
  ics.nscai.syph <- getICsForClusters(read.tree('tree_nscai_excl_20010701_20170501'), 7, 'nscai','syph_after', 'binomial', 'Lab-confirmed syphilis', BIC)
  ics.nscai.num.syph <- getICsForClusters(read.tree('tree_nscai_excl_20010701_20170501'), 7,  'nscai','num_syph_after', 'poisson', 'Number of lab-confirmed syphilis episodes', BIC)
  dfs <- list(ics.nscai.nscai, ics.nscai.std, ics.nscai.syph, ics.nscai.num.std, ics.nscai.num.syph); name <- 'nsCAI'
  }


if (behaviour == 'nsp'){
  ics.nsp.nsp <- getICsForClusters(read.tree('tree_nsp_excl_20010701_20170501'), 10, 'nsp', 'nsp_after', 'binomial', 'nsP', BIC)
  ics.nsp.syph <- getICsForClusters(read.tree('tree_nsp_excl_20010701_20170501'), 10, 'nsp','syph_after', 'binomial', 'Lab-confirmed syphilis', BIC)
  ics.nsp.num.syph <- getICsForClusters(read.tree('tree_nsp_excl_20010701_20170501'), 10,  'nsp','num_syph_after', 'poisson', 'Number of lab-confirmed syphilis episodes', BIC)
  ics.nsp.std <- getICsForClusters(read.tree('tree_nsp_excl_20010701_20170501'), 10, 'nsp','std_after', 'binomial', 'Nurse/physician-reported STIs', BIC)
  ics.nsp.num.std <- getICsForClusters(read.tree('tree_nsp_excl_20010701_20170501'), 10,  'nsp','num_std_after', 'poisson', 'Number of nurse/physician-reported STIs', BIC)
  dfs <- list(ics.nsp.nsp, ics.nsp.std,  ics.nsp.syph, ics.nsp.num.std, ics.nsp.num.syph); name <- 'nsP'
}


# ------ arrange output ---------
df <- do.call(rbind.data.frame, dfs)
df_m <- melt(df, id.vars = c('outcome', 'no_clusters'))
df_m$outcome <- factor(df_m$outcome, levels=c("Nurse/physician-reported STIs", "Lab-confirmed syphilis", name,
                                              "Number of nurse/physician-reported STIs", "Number of lab-confirmed syphilis episodes"))


df_plrt <- df_m[df_m$variable == 'p_LRT',]
df_m <- df_m[df_m$variable != 'p_LRT',]
df_plrt$value <- cut(df_plrt$value, c(0, 0.001, 0.01, 0.05, 0.1, 1), right=FALSE, labels=c('***', '**', '*', '.', ''))
df_plrt$value[is.na(df_plrt$value)] <- ''


# ----------- plot --------------
lim <- 11
mPalette <- economist_pal(fill = TRUE)(9)[c(2, 5, 7, 8)]
t <- theme_minimal() + theme(legend.position = c(0.83, 0.25),
                             legend.text = element_text(size = 7),
                             legend.spacing.x = unit(0.2, 'line'),
                             legend.key.size = unit(0.6, 'line'),
                             legend.title = element_blank(),
                             axis.ticks = element_line(size = 0.25),
                             aspect.ratio = 0.5,
                             panel.grid = element_blank(),
                             plot.title = element_text(hjust = 0.5))

bicplots <- ggplot(df_m[df_m$no_clusters <= lim, ], 
                   aes(x = no_clusters, y = value, colour = variable, fill = variable)) +
  geom_line(size = 1) +  geom_point() +
  facet_wrap(~ outcome, scales = 'free_y') + 
  
  # labels and shaping
  scale_color_manual(values = mPalette, 
                     labels = c('Clusters + age', 
                                paste(name, 'before cut-off + syphilis before cut-off + age'), 
                                paste('Clusters +', name, 'before cut-off + syphilis before cut-off + age'))) +
  scale_fill_manual(values = mPalette, 
                    labels = c('Clusters + age', 
                               paste(name, 'before cut-off + syphilis before cut-off + age'), 
                               paste('Clusters +', name, 'before cut-off + syphilis before cut-off + age'))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.25, 0.05), add=c(0,0))) +
  geom_text(data = df_plrt[df_plrt$no_clusters <= lim,], 
            mapping = aes(x = no_clusters, y = - Inf,  label = value), size = 3, colour = 'black', vjust = -1) +
  guides(fill = FALSE) +
  labs(y = 'BIC', x = 'Number of clusters',
       caption = "\nAsterisks mark p value of LRT of full model without clusters vs. with clusters. *** <0.001; ** <0.01; * <0.05; . <0.1") + t +
  
  # adds "axis lines" to non-border facets
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

bicplots


try(dev.off(), silent = T)
png(paste0(format.Date(today(), '%Y%m%d'), '_stis_noclusters_', name,'.png'), width = 25, height = 10.5, units = 'cm', res = 500)
bicplots
dev.off()

