library(ggsignif)
library(tidyverse)
library(R.matlab)
source('util.R')


#### LOAD FILE USED IN CLUSTERING!!
data_used<-read_csv("data/ADHD_all_features.csv")
results_dir <- "results/ADHD_all_features/"

include_control = FALSE

#### LOAD CLUSTERED DATA
#Open mat file in results directory
clustered_data<- readMat(list.files(results_dir, pattern = "\\.mat$", full.names = TRUE))

#Find number of clusters used
number_clusters <- ncol(clustered_data$CIDX)
#Create columns to append to data
cluster_columns <- paste("cluster", 2:ncol(clustered_data$CIDX), sep = "_")

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


#Append clustered data
for (i in 1:(number_clusters - 1)) {
  data_used[, cluster_columns[i]] <- clustered_data$CIDX[, i + 1]
}


#Combine clinical and non clinical features + remove control patients
data_used <- data_used %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)

#Remove control from analysis
if (include_control == FALSE){
  data_used <- data_used%>%
    filter(dx.original !='CN') 
}


############## RUN HYDRA CLUSTERS ##############


#Run feature analysis on clustered data
feature_analysis(data_used, results_dir)



  
