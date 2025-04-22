library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)

##GLOBAL VS ALL FEATURES
#Directory to save results
saving_dir <- 'results/ADHD_global_all_alluvial.png'

#All features directory
all_features_dir <- 'results/ADHD_all_features'
all_features_clustered <- readMat(list.files(all_features_dir, pattern = "\\.mat$", full.names = TRUE))

#Global features directory
global_features_dir <- 'results/ADHD_global'
global_features_clustered <- readMat(list.files(global_features_dir, pattern = "\\.mat$", full.names = TRUE))

#Check IDs are the same
table(all_features_clustered$ID == global_features_clustered$ID)

#Load data
data_global <- read_csv('data/ADHD_global.csv')
data_all <- read_csv('data/ADHD_all_features.csv')

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


###################CONTROL VS NO CONTROL 
main_dir <-'results/ASD_ADHD_5-10/tsne_medoids/'

control_true <- read_csv(file.path(main_dir, 'control_TRUE', 'tsne_medoids_clustering.csv'))
control_false <- read_csv(file.path(main_dir, 'control_FALSE', 'tsne_medoids_clustering.csv'))

#Find only patient IDs
used_ID <- control_false$ID

control_true_crop <- control_true %>%
  filter(ID %in% used_ID)

combined_data <- data.frame(
  control_TRUE = control_true_crop$cluster_tsne,
  control_FALSE = control_false$cluster_tsne
)

NMI_clustering <-NMI(combined_data$control_TRUE, combined_data$control_FALSE)
ggplot(combined_data,
       aes(y = stat(count), axis1 = control_TRUE, axis2 = control_FALSE)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  scale_x_discrete(limits = c("TRUE", "FALSE"), expand = c(.05, .05)) +
  ggtitle(paste('NMI ', NMI_clustering))

ggsave(file.path(main_dir, 'control_comparison.png'))

print(NMI_clustering)

