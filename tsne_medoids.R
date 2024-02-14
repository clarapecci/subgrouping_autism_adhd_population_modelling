library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
library(ggplot2)
library(orca)

source("util.R")

#Load data from file 
patient_data<-read_csv("data/ASD_global_features.csv")

#Specify results directory 
results_dir <- 'results/tsne_medoids/new_script_asd_global'

#Choose between clustering just patients or patients + control
include_control = TRUE 

######################### SCRIPT #################
#Create results directory if it doesn't exist 
if (!file.exists(results_dir)){
  dir.create(results_dir)
}

#Delete controls if requested
if (include_control ==FALSE){
  patient_data <- patient_data %>%
    filter(group != '-1')
}

#Select features
feature_data <- patient_data %>%
  select(contains('.q.wre')) 

#Carry out clustering
tsne_clustering<- tsne_k_medoids(feature_data[1:100, ])

#Concatenate ID with embeddings and clustering
output_data <- patient_data%>%
  select(ID)%>%
  mutate(embeddings_1 = tsne_clustering[[1]][,1 ])%>%
  mutate(embeddings_2 = tsne_clustering[[1]][,2])%>%
  mutate(cluster_tsne= tsne_clustering[[2]][, 1])

#Write output data
write_csv(output_data, file.path(results_dir, 'tsne_medoids_clustering.csv'))

#Create plot
filename <- file.path(results_dir, 'tsne_medoids_plot.png')
ggplot(output_data, aes (x = embeddings_1, y = embeddings_2, color = clustering)) +
  geom_point()
ggsave(filename)

