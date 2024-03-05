library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
library(ggplot2)
source("util.R")


#Load data from file 
patient_data<-read_csv("data/ASD_global_features.csv")

#Specify directory to save results
results_dir <- 'results/ASD_global/'

#Choose dimensionality reduction method
method ='tsne'

#Specify directory to save clustering results
saving_dir <- file.path(results_dir, paste0(method, '_medoids'))

#Choose whether to run all clusters from k = 2 to 10
run_all = TRUE

#Choose between clustering just patients or patients + control
include_control = FALSE 

#Create new directory for all clustering solutions
if (run_all==TRUE){
  saving_dir= file.path(saving_dir, 'run_all')
  
}

######################### SCRIPT #################
#Create results directory if necessary
if (!file.exists(saving_dir)){
  dir.create(saving_dir)
}

#Delete controls if requested
if (include_control ==FALSE){
  patient_data <- patient_data %>%
    filter(group != '-1')
}

#Create directory for results
saving_dir <- file.path(saving_dir, paste0('control_', include_control))
if (!file.exists(saving_dir)){
  dir.create(saving_dir)
}

#Select features
feature_data <- patient_data %>%
  select(ID, contains('.q.wre')) 

#Carry out clustering
tsne_clustering<- tsne_k_medoids(data = feature_data, method = method, results_dir  = saving_dir, run_all = run_all)



