library(tidyverse)
library(R.matlab)
library(rstatix)
library(perm)
library("dplyr")
source('util.R')

#For the optimal clustering -> calculates significance with respect to control!


#Specify permutation type
perm.opt <- permControl(nmc=10000,seed=123451, setSEED = TRUE)

#### LOAD FILE USED IN CLUSTERING!!
data_used <- read_csv('data/ASD_ADHD_all.csv')
results_dir <- "results/ASD_ADHD_all"

#Directory to save results
main_dir <- file.path(results_dir, 'significance_analysis_HYDRA_perm')

#### LOAD CLUSTERED DATA
#Open mat file in results directory
clustered_data<- readMat(list.files(results_dir, pattern = "\\.mat$", full.names = TRUE))

#Check significant features of optimal cluster!
highest_ARI_cluster <- which.max(clustered_data$ARI)

#Combine all clinical features with data
combined_data <- append_all_features(data_used)

#Combine with clustered data
combined_data <- combined_data %>%
  mutate(cluster_assignment = clustered_data$CIDX[, highest_ARI_cluster])

#Obtain list of features
feature_columns <- combined_data %>%
  select(-c(ID, group), -contains('cluster'))  %>%
  colnames()

#UNCOMMENT RUN TEST ON SMALLER SAMPLE OF DATA
# combined_data <- combined_data%>%
#    select(c(ID, GMVTransformed.q.wre, sGMVTransformed.q.wre, VentriclesTransformed.q.wre, meanCT2Transformed.q.wre, cluster_assignment))

#List of clusters  
clusters <- levels(factor(combined_data$cluster_assignment))

#Remove control group from list of clusters
clusters <- clusters[clusters!='-1']  

#Create empty datafrane to store Cohen's d
cohend_df <- data.frame(matrix( nrow = length(feature_columns), ncol = 4))
colnames(cohend_df)<-c('features', 'control_1', 'control_2', '1_2')
cohend_df[1] <- feature_columns

########START LOOP ###########

#Iterate over clusters and control
for (cluster in clusters){
  
  
  #Create directory to save 
  cluster_dir <- file.path(main_dir, paste0(cluster, '_control'))
  if (!file.exists(cluster_dir)){
    dir.create(cluster_dir)
  }
  
  #Create column for Cohen's d csv
  current_column <-paste0('control_', cluster)
  
  #Iterate over features
  for (i in feature_columns){
    
    #If entire column is NA - skip
    if (all(is.na(combined_data[[i]]))){
      next
    }
    
    #Remove NA in given feature
    data_to_analyse <- combined_data[!is.na(combined_data[[i]]), ]%>%
      filter(cluster_assignment == cluster | cluster_assignment == -1)
    
    
    if (i =='sex'|i =='site' | (i =='dx.original' && nlevels(factor(combined_data$dx.original))>1 )){

      stat_test <- chisq.test(data_to_analyse[[i]], data_to_analyse$cluster_assignment)
      ggplot(data_to_analyse, aes(fill=.data[[i]], x = cluster_assignment)) +
        geom_bar(position="dodge", stat="count") +
        ggtitle(paste('X^2', stat_test$statistic, 'p value ', round(stat_test$p.value, digits = 6)))

    }else{

    #Permutation test
    stat_test <-  permTS(data_to_analyse[[i]][data_to_analyse$cluster_assignment==-1], data_to_analyse[[i]][data_to_analyse$cluster_assignment==cluster], paired=FALSE, alternative="two.sided", method ="exact.mc", control = perm.opt)
    print(paste(i, stat_test$p.value))
    ggplot(data_to_analyse,
           aes(x = factor(cluster_assignment), y = .data[[i]])) +
      xlab('Clusters') +
      geom_boxplot() +
      geom_point(aes(color =.data[[i]]),  position = position_jitter(width = 0.2), alpha = 0.5) +
       scale_color_gradient(low = "blue", high = "red")  +
      ggtitle('p value ', round(stat_test$p.value, digits = 6))
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
    
    #Find Cohen's d
    if (i !='sex'&& i !='site' && i !='dx.original'){
    find_d <- effsize::cohen.d(data_to_analyse[[i]][data_to_analyse$cluster_assignment==-1], data_to_analyse[[i]][data_to_analyse$cluster_assignment==cluster])
    cohend_df[which(cohend_df$features == i), current_column] <- find_d$estimate
    }
    
  }
  
}



#If two clusters -> check against each other
#Select only control and given cluster
no_control <- combined_data %>%
  filter(cluster_assignment != '-1')

#Create directory to save
cluster_dir <- file.path(main_dir, '1_2')
if (!file.exists(cluster_dir)){
  dir.create(cluster_dir)
}

current_column <-paste0('1_2')
for (i in feature_columns){

  #If entire column is NA - skip
  if (all(is.na(combined_data[[i]]))){
    next
  }

  #Remove NA in given feature
  data_to_analyse <- combined_data[!is.na(combined_data[[i]]), ]%>%
    filter(cluster_assignment != -1)

  if (i =='sex'|i =='site' | (i =='dx.original' && nlevels(factor(data_to_analyse$dx.original))>1 )){

    #Calculate chi squared test for non clinical features
    stat_test <- chisq.test(data_to_analyse[[i]], data_to_analyse$cluster_assignment)
    ggplot(data_to_analyse, aes(fill=.data[[i]], x = cluster_assignment)) +
      geom_bar(position="dodge", stat="count") +
      ggtitle(paste('X^2', stat_test$statistic, 'p value ', round(stat_test$p.value, digits = 6)))

  }else if (i !='sex'&& i !='site' && i !='dx.original'){

     #Perform permutation test with clinical features
    stat_test <-  permTS(data_to_analyse[[i]][data_to_analyse$cluster_assignment==1], data_to_analyse[[i]][data_to_analyse$cluster_assignment==2], paired=FALSE, alternative="two.sided", method ="exact.mc", control = perm.opt)
  ggplot(data_to_analyse,
         aes(x = factor(cluster_assignment), y = .data[[i]])) +
    xlab('Clusters') +
    geom_boxplot() +
    geom_point(aes(color =.data[[i]]),  position = position_jitter(width = 0.2), alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red")  +
    ggtitle('p value ', round(stat_test$p.value, digits = 6))

     #Calculate Cohen's d
    find_d <- effsize::cohen.d(data_to_analyse[[i]][data_to_analyse$cluster_assignment==1], data_to_analyse[[i]][data_to_analyse$cluster_assignment==2])
    cohend_df[which(cohend_df$features == i), current_column] <- find_d$estimate
    
  }
  #If p value is significant, include 'significant' directory
  if (!is.na(stat_test$p.value) & stat_test$p.value < 0.005){
    saving_dir <- file.path(cluster_dir, "sig")
  } else{
    saving_dir <- cluster_dir
  }

  # #Attach feature name to saving directory
  filename <- file.path(saving_dir, paste(i, ".png", sep =""))

  #Create plot
  ggsave(filename)



}


#Save Cohen's d csv
write_csv(cohend_df, file.path(main_dir, 'Cohens_d.csv'))




