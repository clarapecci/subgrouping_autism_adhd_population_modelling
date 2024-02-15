library(ggsignif)
library(tidyverse)

source('util.R')


#### LOAD FILE USED IN CLUSTERING!!
data_used<-read_csv("data/ASD_global_features.csv")
results_dir <- "results/tsne_medoids/ASD_all_features/control_TRUE"

#### LOAD CLUSTERED DATA
clustered_data <- read_csv(file.path(results_dir, "tsne_medoids_clustering.csv"))

###### COMBINE FEATURES IN SINGLE DATAFRAME
#Load all clinical data + original features
clinical_data <-read_csv("data/clinical_data_ID.csv")[-c(1)]
all_features_data <- read_csv("data/combat_centile_ID.csv")[-c(1)]

#Check IDs used in clustering to extract desired features from original + clinical data
used_ID <- clustered_data$ID

#Extract correct clinical features
clinical_features <- clinical_data %>% 
  filter(ID %in% used_ID) %>%
  select (-c(participant, site, dx.original, sex, ID))

#Extract non clinical features not used in clustering
non_clinical_features <- all_features_data %>% 
  filter(ID %in% used_ID) %>%
  select(site, age, sex, IQ, dx.original, dx.original)

#Combine clinical and non clinical features
data_used <- data_used %>%
  filter(ID %in% used_ID) %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)



############## RUN FOR TSNE CLUSTERING ########

#Append clusters 
tsne_feature_data <- data_used %>%
  mutate(cluster = clustered_data$cluster_tsne)%>%
  filter(dx.original !='CN')

#Run feature analysis
feature_analysis(tsne_feature_data, results_dir)
