library(tidyverse)
library(filenamer)
library(tidyverse)
library(plotly)
library(Rtsne)
library(fpc)
library(dplyr)
library(cluster)

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

#Function to carry out tSNE dimensionality reduction and clustering with respect to first two embeddings. 
#Input: data frame of only features
#Output: first two embedding, cluster assignment 
tsne_k_medoids <- function(data){
  # run tsne to reduce to 2 dimensions
  sm.tsne <- Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
  # extract the distance matrix from the tsne
  t.dist <- as.matrix(dist(sm.tsne$Y))
  # run partitioning around medoids with silhouette estimation to get the number of optimal clusters
  pamk.best <- pamk(t.dist)
  # run PAM with that number of clusters
  pam.res <- pam(t.dist, pamk.best$nc)
  # put the cluster in a separate variable
  groups <- as.data.frame(pamk.best$pamobject$clustering)
  groups <- as.data.frame(pam.res$clustering)
  
  return(list(sm.tsne$Y, groups))
  
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
  
  
  #Create directory to save figures
  if (!file.exists(results_directory)){ 
    dir.create(results_directory)
  }
  
  #Iterate over number of clusters
  for (x in cluster_list){
    
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
      
      if (i =='sex'|i =='site' | (i =='dx.original' & nlevels(factor(clustered_data_frame$dx.original))>1 )){
        
        stat_test <- chisq.test(clustered_data_frame[[i]], clustered_data_frame[[x]])
        ggplot(clustered_data_frame, aes(fill=.data[[i]], x = .data[[x]])) +
          geom_bar(position="dodge", stat="count") +
          ggtitle(paste('X^2', stat_test$statistic, 'p value ', round(stat_test$p.value, digits = 6)))
        
      }else{

      #Carry out statistical test according to variable type
      stat_test <- kruskal.test(clustered_data_frame[[i]] ~clustered_data_frame[[x]],  data = clustered_data_frame)
      ggplot(clustered_data_frame,
             aes(x = factor(.data[[x]]), y = .data[[i]])) +
        xlab('Clusters') +
        geom_boxplot() +
        ggtitle('p value ', stat_test$p.value)
      
      }
      
      
      #If p value is significant, include 'significant' directory
      if (!is.na(stat_test$p.value) & stat_test$p.value < 0.005){
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
