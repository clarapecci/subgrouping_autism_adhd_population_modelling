library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)


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

#Append clusters to data - no need to include first column as that distinguishes control from patient
combined_data <- data_used%>%
  mutate(cluster_2 = clustered_data$CIDX[,2])%>%
  mutate(cluster_3 = clustered_data$CIDX[,3])%>%
  mutate(cluster_4 = clustered_data$CIDX[,4])%>%
  mutate(cluster_5 = clustered_data$CIDX[,5])%>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)%>%
  filter(dx.original !='CN')

##############ALLUVIAL PLOT ###########

ggplot(combined_data,
       aes(y = stat(count), axis1 = cluster_2, axis2 = cluster_3, axis3 = cluster_4, axis4 = cluster_5)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  scale_x_discrete(limits = c("cluster_2", "cluster_3", "cluster_4", "cluster_5"), expand = c(.05, .05)) 
  #geom_label_repel(stat = "stratum",  box.padding = 0.5, point.padding = 0.1, position = position_nudge_repel(y = 0.5), aes(label = after_stat(stratum)))  



