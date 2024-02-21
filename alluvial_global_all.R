library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)

#Directory to save results
saving_dir <- 'results/ASD_global_all_alluvial.png'

#All features directory
all_features_dir <- 'results/ASD_all_features'
all_features_clustered <- readMat(list.files(all_features_dir, pattern = "\\.mat$", full.names = TRUE))

#Global features directory
global_features_dir <- 'results/ASD_global/'
global_features_clustered <- readMat(list.files(global_features_dir, pattern = "\\.mat$", full.names = TRUE))

#Check IDs are the same
table(all_features_clustered$ID == global_features_clustered$ID)

#Load data
data_global <- read_csv('data/ASD_global_features.csv')
data_all <- read_csv('data/ASD_all_features.csv')

combined_data <- data_all%>%
  mutate(global_clust = global_features_clustered$CIDX[,2] )%>%
  mutate(all_clust = all_features_clustered$CIDX[,2])%>%
  filter(group != -1)
  
NMI_clustering <-NMI(combined_data$global_clust, combined_data$all_clust)
ggplot(combined_data,
       aes(y = stat(count), axis1 = global_clust, axis2 = all_clust)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  scale_x_discrete(limits = c("global", "all"), expand = c(.05, .05)) +
  ggtitle(paste('NMI ', NMI_clustering))
ggsave(saving_dir)
