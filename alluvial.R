library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)
library(aricode)
# TO DO IN SCRIPT:
#   - add sustain


#### LOAD FILE USED IN CLUSTERING!!

#Data used
data_used<-read_csv("data/ASD_global_features.csv")

#Directory where clustering assignments is saved
results_dir <- "results/ASD_global/"

#Specify directory to save new plots
saving_dir <-file.path(results_dir, 'test')

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

#Concatenate with cluster assignment
for (i in 1:(number_clusters - 1)) {
  data_used[, cluster_columns[i]] <- clustered_data$CIDX[, i + 1]
}

#Combine with features
combined_data <- data_used %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)%>%
  filter(dx.original !='CN')
  
##############ALLUVIAL PLOT ###########

filename <- file.path(saving_dir, 'HYDRA_clusters.png')

#Create list of axis based on cluster columns
axis_list <- map(setNames(cluster_columns, paste0("axis", seq_along(cluster_columns))), ~ sym(.x))

# Create the ggplot with axis_list
ggplot(combined_data, aes(y = stat(count), !!!axis_list)) +
      geom_alluvium(width = 1/12) +
       geom_stratum(width = 1/12, fill = "black", color = "grey") +
       scale_x_discrete(limits = unlist(cluster_columns), expand = c(.05, .05))
ggsave(filename)

###### CHECK FOR TSNE + umap
#Check if directory contains tsne file 
if (file.exists(file.path(results_dir, 'tsne_medoids'))){
  
  #Open tsne clustering with no control
  tsne_clustering <- read_csv(file.path(results_dir, 'tsne_medoids', 'control_FALSE', 'tsne_medoids_clustering.csv'))
  
  #Check IDs are the same across HYDRA and tSNE
  stopifnot(tsne_clustering$ID == combined_data$ID)
  
  #Find number of clusters 
  cluster_number_tsne <- nlevels(factor(tsne_clustering$cluster_tsne))
  
  #Find equivalent number of clusters in HYDRA
  cluster_number_tsne <- paste0('cluster_', cluster_number_tsne)
  
  #If number of clusters match
  if (cluster_number_tsne %in% colnames(combined_data)){
    
    #Concatenate tsne-medoids data with HYDRA
    combined_data <- combined_data %>%
      mutate(cluster_tsne = tsne_clustering$cluster_tsne)
    
    #Create directory to save new plot
    compare_plot_dir <- file.path(saving_dir, 'comparison_plot.png')
    
    #Calculate mutual information
    NMI_clustering_tsne <-NMI(combined_data[[cluster_number_tsne]], combined_data$cluster_tsne)
   
    
  #Check for umap clustering  
  if (file.exists(file.path(results_dir, 'umap_medoids', 'control_FALSE'))){
    
    #Open tsne clustering with no control
    umap_clustering <- read_csv(file.path(results_dir, 'umap_medoids', 'control_FALSE', 'umap_medoids_clustering.csv'))
    
    #Check IDs are the same across HYDRA and tSNE
    stopifnot(umap_clustering$ID == combined_data$ID)
    
    #Find number of clusters 
    cluster_number_umap <- nlevels(factor(umap_clustering$cluster_tsne))
    
    #Find equivalent number of clusters in HYDRA
    cluster_number_umap <- paste0('cluster_', cluster_number_umap)
    
    #If number of clusters match
    if (cluster_number_umap %in% colnames(combined_data)){
      
      #Concatenate tsne-medoids data with HYDRA
      combined_data <- combined_data %>%
        mutate(cluster_umap = umap_clustering$cluster_tsne)
      
      #Calculate mutual information
      NMI_clustering_umap <-NMI(combined_data[[cluster_number_tsne]], combined_data$cluster_umap)
      
      ggplot(combined_data,
             aes(y = stat(count), axis1 = .data[[cluster_number_tsne]], axis2 = cluster_tsne, axis3 = cluster_umap)) +
        geom_alluvium(width = 1/12) +
        geom_stratum(width = 1/12, fill = "black", color = "grey") +
        scale_x_discrete(limits = c("HYDRA", "TSNE", "UMAP"), expand = c(.05, .05)) +
        ggtitle(paste('NMI ', NMI_clustering_tsne, NMI_clustering_umap))
      
      ggsave(compare_plot_dir)
    
  }
    }else{
    
    #Alluvial plot across techniques
    ggplot(combined_data,
           aes(y = stat(count), axis1 = .data[[cluster_number_tsne]], axis2 = cluster_tsne)) +
      geom_alluvium(width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      scale_x_discrete(limits = c("HYDRA", "TSNE"), expand = c(.05, .05)) +
      ggtitle('NMI ', NMI_clustering_tsne)
    
    ggsave(compare_plot_dir)
    
  } 
  }
  
  
}
