library(ggsignif)
library(tidyverse)
library(R.matlab)

source('util.R')

#### LOAD FILE USED IN CLUSTERING!!
data_used<-read_csv("data/ASD_all_features.csv")
results_dir <- "results/12-02/ASD_all_features/"

#### LOAD CLUSTERED DATA
clustered_data <- readMat(file.path(results_dir, "ASD_all_features.mat"))

###### COMBINE FEATURES IN SINGLE DATAFRAME
#Load all clinical data + original features
clinical_data <-read_csv("data/clinical_data_ID.csv")[-c(1)]
all_features_data <- read_csv("data/combat_centile_ID.csv")[-c(1)]

#Check IDs used in clustering to extract desired features from original + clinical data
used_ID <- data_used$ID

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
  
#Append clusters to data - no need to include first column as that distinguishes control from patient
combined_data <- data_used%>%
   mutate(cluster_2 = clustered_data$CIDX[,2])%>%
   mutate(cluster_3 = clustered_data$CIDX[,3])%>%
   mutate(cluster_4 = clustered_data$CIDX[,4])%>%
   mutate(cluster_5 = clustered_data$CIDX[,5])%>%
   cbind(clinical_features)%>%
   cbind(non_clinical_features)%>%
   filter(dx.original !='CN')


############## RUN FUNCTION ##############

#Run feature analysis on clustered data
feature_analysis(combined_data, results_dir)


############## RUN FOR TSNE CLUSTERING ########
tsne_clustered <- read_csv('results/tsne_medoids/ASD/ASD_data.csv')%>%
  filter(ID %in% used_ID) 
tsne_feature_data <- data_used %>%
  mutate(cluster = tsne_clustered$tsne_cluster)

tsne_dir <- 'results/tsne_medoids/ASD/'
feature_analysis(tsne_feature_data, tsne_dir)

  
