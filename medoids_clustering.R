library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
library(ggplot2)


source("util.R")

home_dir <- getwd()

#Load data from file 
patient_data<-read_csv(file.path(home_dir, "data/ASD_all_features.csv"))

#Specify directory to save results
results_dir <- file.path(home_dir, 'results/ASD_all_features/')

#Choose dimensionality reduction method
method ='umap'

#Choose between clustering just patients or patients + control
include_control = TRUE 

results_dir <- file.path(results_dir, paste0(method, '_medoids'))
######################### SCRIPT #################
#Create results directory if necessary
if (!file.exists(results_dir)){
  dir.create(results_dir)
}

#Delete controls if requested
if (include_control ==FALSE){
  patient_data <- patient_data %>%
    filter(group != '-1')
}

#Create directory for results
results_dir <- file.path(results_dir, paste0('control_', include_control))
if (!file.exists(results_dir)){
  dir.create(results_dir)
}

#Select features
feature_data <- patient_data %>%
  select(contains('.q.wre')) 

#Carry out clustering
tsne_clustering<- tsne_k_medoids(data = feature_data, method = method)

#Concatenate ID with embeddings and clustering
output_data <- patient_data%>%
  select(ID)%>%
  mutate(embeddings_1 = tsne_clustering[[1]][,1 ])%>%
  mutate(embeddings_2 = tsne_clustering[[1]][,2])%>%
  mutate(cluster_tsne= tsne_clustering[[2]][, 1])

#Write output data
write_csv(output_data, file.path(results_dir, paste0(method, '_medoids_clustering.csv')))

#Create plot
filename <- file.path(results_dir, paste0(method, '_medoids_plot.png'))
ggplot(output_data, aes (x = embeddings_1, y = embeddings_2, color = cluster_tsne)) +
  geom_point()
ggsave(filename)

