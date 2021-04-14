# behavioural-clusters
Codes used to produce the results in "Unsupervised machine learning identifies long-term behavioural patterns and  predicts future sexual behaviour and sexually transmitted infections among HIV-positive men who have sex with men"

The codes are structured as pictured:
![Code Structure](/code_structure.png)


The files do the following:

Data generation
- __datagen.R__: defines _dataGen_, _dataPrepSHCS_ and _dataPrepZPHI_ functions that create a synthetic data set with variables required for analyses


Clustering 
- __0_functions.R__: defines clustering function _getBehaviouralClusters_, and _plotTrends_ function for plotting nsCAI/nsP trends in clusters, along with a few other useful functions


Analysis
- __analyse_stis_regression.R__: defines _runRegression_ function with regression models to test whether clusters are predictive of nsCAI/nsP and STIs 
- __analyse_stis_roc.R__: defines _getROCValues_ function which generates receiver operator characteristic (ROC) values to test whether clusters are predictive of nsCAI/nsP and STIs 
- __analyse_stis_noclusters__: defines _getICsForClusters_ function which runs regression models for different numbers of behavioural clusters and compares their predictivity for nsCAI/nsP and STIs using information criteria (ICs)
- __analyse_nopartners__: defines _runRegression.nopartners_ functions with mixed effect regression models to test whether clusters are associated with participants' number of partners
- __analyse_nopartners_noclusters__: defines _getICsForClusters.nopartners_ function which runs regression models for different numbers of behavioural clusters to and compares their association with participants' number of partners using information criteria (ICs)


Plotting
- __plot_stis_forest.R__: calls _runRegression_ from analyse_stis_regression.R, creates forest plots
- __plot_stis_bic.R__: calls _runRegression_ from analyse_stis_regression.R, creates bar plot with BIC values
- __plot_stis_roc.R__: calls _getROCValues_ from analyse_stis_roc.R, creates ROC curves and corresponding AUC bar plots
- __plot_stis_noclusters.R__: calls _getICSForClusters_ from analyse_stis_noclusters.R, creates line plots with BIC values
- __plot_nopartners_forest.R__: calls _runRegression.nopartners_ from analyse_nopartners.R, creates forest plots
- __plot_nopartners_noclusters.R__: calls _getICsForClusters.nopartners_ from analyse_nopartners_noclusters.R, creates line plots wth BIC values
- __table_stis_ics__: calls _runRegression_ from analyse_stis_regression.R, creates table with p_LRTs and information criteria
