library(umap)
library(tidyverse)
library(plotly)
library(dplyr)
library(tidytext)
source("util.R")

#Load data
##-- directly from file
patient_data<-read_csv("data/ASD_all_features.csv")

#Select features
feature_data <-patient_data %>%
  select(contains('q.wre'))

#Carry out UMAP
umap_data <- umap(feature_data)
embeddings <- umap_data[["layout"]] 
embeddings <- data.frame(embeddings)

#Load all clinical data + original features
clinical_data <-read_csv("data/clinical_data_ID.csv")[-c(1)]
all_features_data <- read_csv("data/combat_centile_ID.csv")[-c(1)]

#Check IDs used in clustering to extract desired features from original + clinical data
used_ID <- patient_data$ID

#Extract correct clinical features
clinical_features <- clinical_data %>% 
  filter(ID %in% used_ID) %>%
  select (-c(participant, site, dx.original, sex, ID))

#Extract non clinical features 
non_clinical_features <- all_features_data %>% 
  filter(ID %in% used_ID) %>%
  select(site, age, sex, IQ, dx.original, dx.original)

#Concatenate UMAP data with original data for plotting
data_plus_umap  <-  patient_data %>%
  mutate(X1 = embeddings$X1) %>%
  mutate(X2 = embeddings$X2)

#Concatenate with features
data_plus_umap <- data_plus_umap %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)

#Change color variable with any feature of interest!  
plot_ly(data_plus_umap, x = ~X1, y = ~X2, color = ~age, colors = c('#636EFA', '#EF553B', '#00CC96'), type = 'scatter', mode = 'markers')
