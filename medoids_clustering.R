library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
library(ggplot2)
source("util.R")


#Load data from file 
data_to_cluster_dir <- "/your/path/here/"  # <-- EDIT THIS PATH
data_to_cluster<-read_csv(data_to_cluster_dir)

#Specify directory to save results
results_dir <- "/your/path/here/results"  # <-- EDIT THIS PATH

#Choose whether to save all clusters from k = 2 to 10 (TRUE) or just the optimal (FALSE)
run_all = FALSE

#Choose between clustering just patients or patients + control
include_control = FALSE 

#Choose dimensionality reduction method ('umap' or 'tsne')
method ='umap'

#Specify directory to save clustering results
saving_dir <- file.path(results_dir, paste0(method, '_medoids'))



######################### SCRIPT #################
#Create new directory for all clustering solutions
if (run_all==TRUE){
  saving_dir= file.path(saving_dir, 'run_all')
  
}

#Delete controls if requested
if (include_control ==FALSE){
  data_to_cluster <- data_to_cluster %>%
    filter(group != '-1')
}

#Create directory for results
saving_dir <- file.path(saving_dir, paste0('control_', include_control))

#Data to be clustered should only contain the IDs and the features to be used (in our case, the feature columns we wanted to use all contained '.q.wre' in the name)
feature_data <- data_to_cluster %>%
  select(ID, contains('.q.wre')) 

#Create results directory
if (!file.exists(saving_dir)){
  dir.create(saving_dir, recursive = TRUE)
}

#Carry out clustering
clustering_solution<- k_medoids_clustering(data = feature_data, method = method, results_dir = saving_dir, run_all = run_all)




