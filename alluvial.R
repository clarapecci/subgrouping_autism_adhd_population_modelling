library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)
library(aricode)
library(ggfittext)
source('util.R')


#Load .csv used in clustering
data_used_dir <-  '/your/path/here/' # <-- EDIT THIS PATH
data_used <- read_csv(data_used_dir)

#Specify file where clustered data is saved
results_dir <- "/your/path/here/"  # <-- EDIT THIS PATH

#Specify directory to save new plots
saving_dir <-file.path(results_dir, 'comparisons_plots')

#Create directory to save results
if (!file.exists(saving_dir)){
  dir.create(saving_dir)
}

#### LOAD CLUSTERED DATA
#Open mat file in results directory
clustered_data<- readMat(list.files(results_dir, pattern = "\\.mat$", full.names = TRUE))

#Find number of clusters and create new column names
number_clusters <- ncol(clustered_data$CIDX)
cluster_columns <- paste("cluster", 2:ncol(clustered_data$CIDX), sep = "_")

###### COMBINE FEATURES IN SINGLE DATAFRAME
#Combine all clinical features with data
combined_data <- append_all_features(data_used)

#Find optimal cluster 
highest_ARI_cluster <- which.max(clustered_data$ARI)

#Append all clusters to data
for (i in 1:(number_clusters - 1)) {
  combined_data[, cluster_columns[i]] <- clustered_data$CIDX[, i + 1]
}


combined_data <- combined_data%>%
  filter(dx.original !='CN')
  

  


########### find NMI between solutions in same technique ########### 

#Create list of clusters used
cluster_list <- seq(from = 2, to =dim(clustered_data$CIDX)[2] )

#Create all combinations between clusters
cluster_combs <- combn(cluster_list, 2)

#Empty array to store values
loop_NMI <- numeric(length = dim(cluster_combs)[2])

#Iterate over all cluster combinations
for (i in 1:dim(cluster_combs)[2]){

  first_cluster <- paste0('cluster_', cluster_combs[1, i])
  second_cluster <- paste0('cluster_', cluster_combs[2, i])

  #Calculate NMI between cluster pairs
  loop_NMI[i] <- NMI(combined_data[[first_cluster]], combined_data[[second_cluster]])
}

print(paste('mean', mean(loop_NMI)))
print(paste('sd', sd(loop_NMI)))

#Plot alluvial across clusters in the same method
filename <- file.path(saving_dir, 'HYDRA_clusters.png')

#Create list of axis based on cluster columns
axis_list <- map(setNames(cluster_columns, paste0("axis", seq_along(cluster_columns))), ~ sym(.x))

# Create the ggplot with axis_list
 ggplot(combined_data, aes(y = stat(count),!!!axis_list)) +
       geom_alluvium(width = 1/12) +
        geom_stratum(width = 1/12, fill = "black", color = "grey", spacing = 1000) +
   theme_bw()+
        scale_x_discrete(limits = unlist(cluster_columns), expand = c(.05, .05))
 ggsave(filename, create.dir = TRUE)




###### COMPARE HYDRA TO OTHER METHODS
#Find highest ARI cluster in HYDRA to compare to others
highest_ARI_cluster <-  paste0('cluster_', which.max(clustered_data$ARI) ) 


#Check for umap clustering
if (file.exists(file.path(results_dir, 'umap_medoids', 'run_all', 'control_FALSE', '1'))){

  #Open tsne clustering with no control
  umap_clustering <- read_csv(file.path(results_dir, 'umap_medoids', 'run_all','control_FALSE', '1', 'umap_medoids_clustering.csv'))

  #Check IDs are the same across HYDRA and K-medoids
  stopifnot(umap_clustering$ID == combined_data$ID)

  #Concatenate k-medoids solution with HYDRA
  combined_data <- combined_data %>%
    mutate(cluster_umap = umap_clustering$cluster_umap)
  
  my_colors <- c("red", "blue", "green")
  
  #Calculate mutual information
  NMI_clustering_umap <-round(NMI(combined_data[[highest_ARI_cluster]], combined_data$cluster_umap), 4)
  
  compare_plot_dir <- file.path(saving_dir, 'comparison_plot.png')

  ggplot(combined_data,
         aes(y = stat(count), axis1 = .data[[highest_ARI_cluster]], axis2 = cluster_umap, fill = as.factor(.data[[highest_ARI_cluster]]))) +
    geom_alluvium( width = 1/12 ) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    scale_x_discrete(limits = c("HYDRA", "UMAP"), expand = c(.05, .05)) +
    scale_fill_manual(values = my_colors)  + 
    labs(fill = "HYDRA subgroups") +
    ggtitle(paste('NMI ', NMI_clustering_umap))

  ggsave(compare_plot_dir, create.dir = TRUE)

  }


