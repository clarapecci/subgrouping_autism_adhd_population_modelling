library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
library(ggplot2)
source("util.R")

#Load data from file 
patient_data<-read_csv("data/ASD_ADHD_global_5-10.csv")

#Specify directory to save results
results_dir <- 'results/ASD_ADHD_5-10'

#Choose dimensionality reduction method
method ='umap'

#Choose between clustering just patients or patients + control
include_control = FALSE 

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

#Create new column name depending on method 

new_column <- paste0('cluster_', method)
#Concatenate ID with embeddings and clustering
output_data <- patient_data%>%
  select(ID)%>%
  mutate(embeddings_1 = tsne_clustering[[1]][,1 ])%>%
  mutate(embeddings_2 = tsne_clustering[[1]][,2])%>%
  mutate(!!new_column := tsne_clustering[[2]][, 1])

#Write output data
write_csv(output_data, file.path(results_dir, paste0(method, '_medoids_clustering.csv')))

#Create plot
filename <- file.path(results_dir, paste0(method, '_medoids_plot.png'))
ggplot(output_data, aes (x = embeddings_1, y = embeddings_2, color = !!sym(new_column))) +
  geom_point() + scale_color_gradient(low = "blue", high = "red") 

ggsave(filename)

