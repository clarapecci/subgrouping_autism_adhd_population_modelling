library(tidyverse)
library(filenamer)
library(tidyverse)
library(plotly)
library(Rtsne)
library(fpc)
library(dplyr)
library(cluster)
library(umap)

################## START OF FUNCTION ###################
#Functions takes name of csv file and returns processed data frame where it can select what type of patient (ASD, ADHD, both), in what window age, and what number of samples (crop)
preprocess_data <- function (csv_file, patient, min_age = NULL, max_age = NULL , crop = NULL, saving_dir = NULL ){
  
  #Load original data
  original_data <- read_csv(csv_file)%>%
    select(c(-1))
  
  #Give each subject a numeric ID
  ID_data <-original_data %>%
    mutate(ID = row_number())
  
  #Remove low quality samples
  QC_data<- ID_data %>%
    filter(FSQC <2.5)%>%
    filter(site != 'UCSD')%>% 
    filter(site!='UCSDnew') #Remove UCSD data
  
  
  #Change column name (avoid '.')
  colnames(QC_data)[colnames(QC_data) == 'dx.original'] <- 'diagnostic'  

  
  #Select patient type
  if (patient == 'ASD'){
    patient_data<- QC_data %>%
      filter(diagnostic!='ADHD')%>%
      mutate(group = case_when(diagnostic =='ASD'~ 1,
                               diagnostic =='CN'~ -1,))
  } else if (patient == 'ADHD') {
    patient_data<- QC_data %>%
      filter(diagnostic!='ASD')%>%
      mutate(group = case_when(diagnostic =='ADHD'~ 1,
                               diagnostic =='CN'~ -1,))
  } else if (patient =='both'){
    patient_data<- QC_data %>%
      mutate(group = case_when(diagnostic =='ADHD'~ 1,
                               diagnostic == 'ASD' ~1, 
                               diagnostic =='CN'~ -1,))
  }
  
  #Age window
  #Minimum
  if(!is.null(min_age)){
    patient_data <- patient_data%>%
      filter(age >min_age)
  }
  
  #Maximum
  if(!is.null(max_age)){
    patient_data <- patient_data%>%
      filter(age < max_age)
  }
  
  
  #Crop data
  if (!is.null(crop)){
    patient_data <- patient_data[1:crop, ]
  }
  
  #Save data
  if (!is.null(saving_dir)){
    write.csv(patient_data, file.path(new_data_filename+'covariate.csv'))
  }
  
  return (patient_data)
  
}
############## END OF FUNCTION #####################


################## START OF FUNCTION ###################

#Function to carry out tSNE or UMAP dimensionality reduction and clustering with respect to first two embeddings. 
#Input: data frame of only ID + features, method = 'tsne' or 'umap' to decide on dimensionality reduction type
#Output: first two embedding, cluster assignment 
#If run_all is TRUE -> clustering is run for all clusters k = 1 to 10 // otherwise only output results for ideal clustering
k_medoids_clustering <- function(data, method, results_dir, run_all =FALSE){
  
  #Name of column with optimal clustering
  optimal_column <- paste0('cluster_', method)
  
  #Select features from data
  data<- feature_data %>%
    select(contains('.q.wre'))
  
  #Select the ID from the input data
  clustering_data_frame <- feature_data%>%
    select(ID)
  
  #Apply dimensionality reduction
  if (method =='tsne'){
    print('tSNE method chosen... applying dimensionality reduction')
    # run tsne to reduce to 2 dimensions
    sm.tsne <- Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
    embeddings <- sm.tsne$Y
  } else if (method =='umap'){
    print('UMAP method chosen... applying dimensionality reduction')
    umap_data <- umap(data)
    embeddings <- umap_data$layout
  }else{
    print('Method must be either tsne or umap')
    return ()
  }

  #Add embeddings to output data frame
  clustering_data_frame <- clustering_data_frame%>%
    mutate(embeddings_1 = embeddings[,1 ])%>%
    mutate(embeddings_2 = embeddings[,2])
  
  # extract the distance matrix from the dimensionality reduction
  t.dist <- as.matrix(dist(embeddings))

  #Run for all clustering solutions
  if (run_all ==TRUE){
    
    #Create empty vector to store all Silhouette scores
    criterion <-numeric(10)
    
    #Run clustering for all clusters from k = 1 to k + 10
    for (i in 2:10){
      
      #Find pam object 
      pam_value<- pam(t.dist, i)
      
      print(paste('cluster ', i, 'silhoutte score ', pam_value$silinfo$avg.width))
      current_column <- paste0('cluster_', i)
      
      #Add clustering to
      clustering_data_frame <- clustering_data_frame%>%
        mutate(!!current_column := pam_value$clustering)
      
      #Append Silhouette score to criterion vector
      criterion[i] <- pam_value$silinfo$avg.width
     
      }
    
    #Find clustering with highest silhouette score 
    best_cluster <- paste0('cluster_', which.max(criterion))
    print(paste('Best cluster number: ', best_cluster))
    clustering_data_frame <- clustering_data_frame%>%
      mutate(!!optimal_column := clustering_data_frame[[best_cluster]])
  
  
  }
  
  if (run_all ==FALSE){
    #Only run the clustering for the highest silhouette score 
    
    # run partitioning around medoids with silhouette estimation to get the number of optimal clusters
    print('Finding optimal number of clusters ....')
    pamk.best <- pamk(t.dist, critout = TRUE)
    
    criterion <- pamk.best$crit
    # run PAM with that number of clusters
    print('Running partitioning around medoids...')
    
    #Select optimum clustering
    pam.res <- pamk.best$pamobject
    
    #Concatenate ID with embeddings and clustering
    clustering_data_frame <- clustering_data_frame%>%
      mutate(!!optimal_column := pam.res$clustering)
    }
  
  #Save criterion params
  if (!is.null(results_dir)){
  
  print(paste('Saving results to ', results_dir))
  criterion_save <- file.path(results_dir, 'criterion.csv')
  write.csv(criterion, criterion_save, row.names=FALSE)
  
  #Write output data
  write_csv(clustering_data_frame, file.path(results_dir, paste0(method, '_medoids_clustering.csv')))
  
  #Create plot of embeddings coloured by clustering
  filename <- file.path(results_dir, paste0(method, '_medoids_plot.png'))
  ggplot(clustering_data_frame, aes (x = embeddings_1, y = embeddings_2, color = !!sym(optimal_column))) +
    geom_point() + scale_color_gradient(low = "blue", high = "red") 
  
  ggsave(filename)
  }
}

############## END OF FUNCTION #####################



################## START OF FUNCTION ###################
#Function that takes data frame including clustered data, and makes boxplots for each feature separated into the different clusters.
#If the difference between the clusters is significant, the plots are saved in a 'sig' folder
#Input: data frame including features to be selected and clustered information, directory to save plots in 
feature_analysis <- function(clustered_data_frame, results_directory){
  
  #Find total number of clusters
  cluster_list <- colnames(select(clustered_data_frame, contains('cluster')))
  
  #List of features to iterate over
  feature_columns <- clustered_data_frame %>%
    select(-c(ID, group), -contains('cluster'))  %>%
    colnames()
  
  print(cluster_list)
  #Create directory to save figures
  if (!file.exists(results_directory)){ 
    dir.create(results_directory)
  }
  
  #Iterate over number of clusters
  for (x in cluster_list){
    print(x)
    #Create file directory for given cluster
    cluster_dir <- file.path(results_directory, x)
    if (!file.exists(cluster_dir)){
      dir.create(cluster_dir)
    }
    
    for (i in feature_columns){

      print(i)

      #If entire column is NA - skip
      if (all(is.na(clustered_data_frame[[i]]))){
        next
      }
      
      if (i =='sex'|i =='site' | (i =='dx.original' && nlevels(factor(clustered_data_frame$dx.original))>1 )){
        
        stat_test <- chisq.test(clustered_data_frame[[i]], clustered_data_frame[[x]])
        ggplot(clustered_data_frame, aes(fill=.data[[i]], x = .data[[x]])) +
          geom_bar(position="dodge", stat="count") +
          ggtitle(paste('X^2', stat_test$statistic, 'p value ', round(stat_test$p.value, digits = 6)))
        
      }else if (i !='sex'&& i !='site' && i !='dx.original'){

      #Carry out statistical test according to variable type
      stat_test <- kruskal.test(clustered_data_frame[[i]] ~clustered_data_frame[[x]],  data = clustered_data_frame)
      ggplot(clustered_data_frame,
             aes(x = factor(.data[[x]]), y = .data[[i]])) +
        xlab('Clusters') +
        geom_boxplot() +
        geom_point(aes(color = .data[[i]]), position = 'jitter') +  scale_color_gradient(low = "blue", high = "red")  + 
        ggtitle('p value ', stat_test$p.value)
      
      }
      
      
      #If p value is significant, include 'significant' directory
      if (!is.na(stat_test$p.value) & stat_test$p.value < 0.05){
        saving_dir <- file.path(cluster_dir, "sig")
      } else{
        saving_dir <- cluster_dir
      }

      #Attach feature name to saving directory
      filename <- file.path(saving_dir, paste(i, ".png", sep =""))

      #Create plot
      ggsave(filename)

    }
    
  }
  
}

############## END OF FUNCTION #####################

##############START OF FUNCTION ###################
#Find number and type of significant features for given clustering
find_sig_features <- function(results_dir, data_used, clinical_features, non_clinical_features){
  
  #Check if any features are significant
  sig_dir <- file.path(results_dir, 'sig')
  if (file.exists(sig_dir)){
    
    #Find significant features
    all_files <- list.files(sig_dir, pattern = ".png", full.names = FALSE)
    
    #Obtain name of significant features
    variable_names <- gsub(".png", "", basename(all_files))
    
    #Find anatomical features
    anatomical_features <- data_used%>% select(contains('q.wre')) %>%colnames
    
    #Count anatomical features
    anatomical_count <- length(intersect(variable_names, anatomical_features))
    
    #Count clinical features
    clinical_count <- length(intersect(variable_names, colnames(clinical_features)))
    
    #Count non clinical features
    non_clinical_count <- length(intersect(variable_names, colnames(non_clinical_features)))
    
    #Print number of features 
    print(paste('Significant features in ...', sig_dir))
    print(paste("Anatomical:", anatomical_count))
    print(paste('Clinical:', clinical_count))
    print(paste('Non clinical: ', non_clinical_count))
  }else{
    print(paste('No significant features in ...', sig_dir))
  }
  
  
}


###########START OF FUNCTION ######################
#Concatenate features used in clustering with all clinical and non clinical features for analysis
append_all_features <- function(data_used){
  
  #Load all clinical data
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
    select(site, age, sex, IQ, dx.original, FSQC)
  
  #Append features
  combined_data <- data_used %>%
    cbind(clinical_features)%>%
    cbind(non_clinical_features)
  
  return(combined_data)
}
###########END OF FUNCTION ######################

