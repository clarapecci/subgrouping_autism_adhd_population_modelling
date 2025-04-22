library(ggsignif)
library(tidyverse)
library(R.matlab)
library(ggalluvial)
library(ggrepel)
library(aricode)
library(ggfittext)
# TO DO IN SCRIPT:

#### LOAD FILE USED IN CLUSTERING!!

#Data used
data_used<-read_csv("data/ASD_ADHD_global.csv")

#Directory where clustering assignments is saved
results_dir <- "results/sept_results/ASD_ADHD_global_sept/"

# #Specify directory to save new plots
# saving_dir <-file.path(results_dir, 'comparisons_new_plots')
# 
# #Create directory to save results
# if (!file.exists(saving_dir)){ 
#   dir.create(saving_dir)
# }

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
  select (-c(site, dx.original, sex, ID))

#Extract list of participants for sustain
participants_used <- clinical_features$participant

#Extract non clinical features not used in clustering
non_clinical_features <- all_features_data %>% 
  filter(ID %in% used_ID) %>%
  select(site, age, sex, IQ, dx.original, dx.original)

# #Concatenate with cluster assignment
# for (i in 1:(number_clusters - 1)) {
#   data_used[, cluster_columns[i]] <- clustered_data$CIDX[, i + 1]
# }

highest_ARI_cluster <- which.max(clustered_data$ARI)
#Combine with features
combined_data <- data_used %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)%>%
  mutate(cluster_hydra = clustered_data$CIDX[, highest_ARI_cluster])%>%
  filter(dx.original !='CN')
  

  


########### find NMI between clusters in same technique ########### 
# saving_df <- read_csv('results/NMI_in_technique.csv')
# #Create list of clusters used
# cluster_list <- seq(from = 2, to =dim(clustered_data$CIDX)[2] )
# 
# #Create all combinations between clusters
# cluster_combs <- combn(cluster_list, 2)
# 
# #Empty array to store values 
# loop_NMI <- numeric(length = dim(cluster_combs)[2])
# 
# #Iterate over all cluster combinations
# for (i in 1:dim(cluster_combs)[2]){
#   
#   first_cluster <- paste0('cluster_', cluster_combs[1, i])
#   second_cluster <- paste0('cluster_', cluster_combs[2, i])
#   
#   #Calculate NMI between cluster pairs
#   loop_NMI[i] <- NMI(combined_data[[first_cluster]], combined_data[[second_cluster]])
# }
# 
# print(paste('mean', mean(loop_NMI)))
# print(paste('sd', sd(loop_NMI)))
# 
# #Create new row to add to data frame
# new_row <- tibble('File' =results_dir,
#                       'All_NMI' = list(array(loop_NMI)), 
#                       'Mean' = mean(loop_NMI),
#                       'Std' = sd(loop_NMI))
# 
# #Save new data frame
# updated_df <- rbind(saving_df, new_row)
# write_csv(updated_df, 'results/NMI_in_technique.csv')


############## ALLUVIAL PLOT ###########
#filename <- file.path(saving_dir, 'HYDRA_clusters.png')

#Create list of axis based on cluster columns
#axis_list <- map(setNames(cluster_columns, paste0("axis", seq_along(cluster_columns))), ~ sym(.x))

# # Create the ggplot with axis_list
#  ggplot(combined_data, aes(y = stat(count),!!!axis_list)) +
#        geom_alluvium(width = 1/12) +
#         geom_stratum(width = 1/12, fill = "black", color = "grey", spacing = 1000) +
#    theme_bw()+
#         scale_x_discrete(limits = unlist(cluster_columns), expand = c(.05, .05))
#  #ggsave(filename)
# 



###### CHECK FOR TSNE + umap 
#Check if directory contains tsne file
#Find highest ARI cluster in HYDRA to compare to others
highest_ARI_cluster <- which.max(clustered_data$ARI)

 
# if (file.exists(file.path(results_dir, 'tsne_medoids'))){
# 
#   #Open tsne clustering with no control
#   tsne_clustering <- read_csv(file.path(results_dir, 'tsne_medoids', 'run_all', 'control_FALSE', 'tsne_medoids_clustering.csv'))
# 
#   #Check IDs are the same across HYDRA and tSNE
#   stopifnot(tsne_clustering$ID == combined_data$ID)
# 
#   #Concatenate tsne-medoids data with HYDRA
#   combined_data <- combined_data %>%
#     mutate(cluster_tsne = tsne_clustering$cluster_tsne)
# 
#   #Create directory to save new plot
#   compare_plot_dir <- file.path(saving_dir, 'comparison_plot.png')
# 
#   #Calculate mutual information
#   NMI_clustering_tsne <-round(NMI(combined_data[[highest_ARI_cluster]], combined_data$cluster_tsne), 4)


  #Check for umap clustering
  if (file.exists(file.path(results_dir, 'umap_medoids', 'run_all', 'control_FALSE', '1'))){

    #Open tsne clustering with no control
    umap_clustering <- read_csv(file.path(results_dir, 'umap_medoids', 'run_all','control_FALSE', '1', 'umap_medoids_clustering.csv'))

    #Check IDs are the same across HYDRA and tSNE
    stopifnot(umap_clustering$ID == combined_data$ID)

    #Concatenate tsne-medoids data with HYDRA
    combined_data <- combined_data %>%
      mutate(cluster_umap = umap_clustering$cluster_umap)# %>%
     #  mutate(cluster_umap = case_when(umap_clustering$cluster_umap =='1'~ 2,
     #                                  umap_clustering$cluster_umap =='2'~ 1))%>%
     # mutate(cluster_tsne = case_when(tsne_clustering$cluster_tsne =='1'~ 2,
     #                                tsne_clustering$cluster_tsne =='2'~ 1))
    
    my_colors <- c("red", "blue", "green")
    #Calculate mutual information
    NMI_clustering_umap <-round(NMI(combined_data[[highest_ARI_cluster]], combined_data$cluster_umap), 4)
    #NMI_tsne_umap <- round(NMI(combined_data$cluster_tsne, combined_data$cluster_umap), 4)
    ggplot(combined_data,
           aes(y = stat(count), axis1 = .data[[highest_ARI_cluster]], axis2 = cluster_umap, fill = as.factor(.data[[highest_ARI_cluster]]))) +
      geom_alluvium( width = 1/12 ) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      scale_x_discrete(limits = c("HYDRA", "UMAP"), expand = c(.05, .05)) +
      scale_fill_manual(values = my_colors)  + 
      labs(fill = "HYDRA subgroups") +
      ggtitle(paste('NMI ', NMI_clustering_umap))

    #ggsave(compare_plot_dir)

    }else{

    #Alluvial plot across techniques
    ggplot(combined_data,
           aes(y = stat(count), axis1 = .data[[highest_ARI_cluster]], axis2 = cluster_tsne)) +
      geom_alluvium(width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      scale_x_discrete(limits = c("HYDRA", "TSNE"), expand = c(.05, .05)) +
      ggtitle('NMI ', NMI_clustering_tsne)

    #ggsave(compare_plot_dir)

  }




#Check number of significant features shared between clusters 
# hydra_dir <- list.files(file.path(results_dir, highest_ARI_cluster, 'sig'))
# tsne_dir <- list.files(file.path(results_dir, 'tsne_medoids', 'control_FALSE', 'cluster', 'sig'))
# umap_dir <- list.files(file.path(results_dir, 'umap_medoids', 'control_FALSE', 'cluster', 'sig'))
# 
# print(paste('HYDRA and tSNE: ', length(intersect(hydra_dir, tsne_dir))))
# print(paste('HYDRA and UMAP: ', length(intersect(hydra_dir, umap_dir))))
# print(paste('UMAP and tsne: ', length(intersect(umap_dir, tsne_dir))))
# 


# if (file.exists(file.path(results_dir, 'ADHD_sustain_subtypes_combat_3subgroups.csv'))){
#   
#   #Open tsne clustering with no control
#   sustain_clustering <- read_csv(file.path(results_dir, 'ADHD_sustain_subtypes_combat_3subgroups.csv'))
#   
#   sustain_clustering <- sustain_clustering%>%
#     filter(participant %in% participants_used)%>%
#     filter(Diagnosis!='CN')
#   
#   #Concatenate tsne-medoids data with HYDRA
#   combined_data <- combined_data %>%
#     mutate(cluster_sustain = sustain_clustering$ml_subtype)
#   
#   #Calculate mutual information
#   NMI_clustering_sustain <-NMI(combined_data[[highest_ARI_cluster]], combined_data$cluster_sustain)
#   ggplot(combined_data,
#          aes(y = stat(count), axis1 = cluster_sustain, axis2 = .data[[highest_ARI_cluster]], axis3 = cluster_tsne, axis4 = cluster_umap) )+
#     geom_alluvium(width = 1/12) +
#     geom_stratum(width = 1/12, fill = "black", color = "grey") +
#     scale_x_discrete(limits = c("SUSTAIN", "HYDRA", "TSNE", "UMAP"), expand = c(.05, .05)) +
#     ggtitle(paste('NMI ', NMI_clustering_tsne, NMI_clustering_umap, NMI_tsne_umap, NMI_clustering_sustain))
#   
#   #Create directory to save new plot
#   compare_plot_dir_2 <- file.path(saving_dir, 'comparison_plot_sustain.png')
#   #ggsave(compare_plot_dir_2)
#   
# }


#Plot with all the global features coloured in by HYDRA subgroup
ggplot(combined_data,
       aes(y = stat(count), axis1 =cluster_hydra, axis2 = GMVTransformed.q.wre, axis3 = totalSA2Transformed.q.wre, axis4 = meanCT2Transformed.q.wre, axis5 = sGMVTransformed.q.wre, axis6 = VentriclesTransformed.q.wre, axis7 = WMVTransformed.q.wre, fill = as.factor(cluster_hydra))) +
  geom_alluvium( width = 1/12 ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  scale_x_discrete(limits = c("HYDRA", "GMV", 'SA', 'CT', 'sGMV', 'Ven', 'WMV'), expand = c(.05, .05)) +
  #scale_fill_gradient(low ='blue', high = 'red')
  #scale_fill_manual(values = my_colors)  + 
  labs(fill = "HYDRA subgroups") 
  #ggtitle(paste('NMI ', NMI_clustering_umap))


#Plot with respect to GMVTransformed - to figure out the order of the centiles
ggplot(combined_data,
       aes(y = stat(count), axis2 = GMVTransformed.q.wre, axis4 = meanCT2Transformed.q.wre, fill = GMVTransformed.q.wre)) +
  geom_alluvium( width = 1/12 ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  #scale_x_discrete(limits = c( "GMV", 'SA', 'CT'), expand = c(.05, .05)) +
  #scale_fill_manual(values = my_colors)  +
  #scale_fill_gradient(low = "blue", high = "red")
  #labs(fill = "HYDRA subgroups")
#ggtitle(paste('NMI ', NMI_clustering_umap))


data_to_pivot <- combined_data%>%
  select(c(ID, GMVTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre, cluster_hydra ))
long_data <- data_to_pivot %>%
  pivot_longer(
    cols =  c(GMVTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre), 
    names_to = "Measure", 
    values_to = "Value"
  )

long_data$cluster_hydra <- as.factor(long_data$cluster_hydra)
long_data$Measure <- as.factor(long_data$Measure)
ggplot(long_data,
       aes(x = Measure, alluvium = cluster_hydra,
           fill = cluster_hydra, label = cluster_hydra)) +  
  geom_alluvium(aes(weight = Value), width = 1/12) +  # Use Value to scale flows
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  labs(fill = "HYDRA subgroups") +
  theme_minimal()
