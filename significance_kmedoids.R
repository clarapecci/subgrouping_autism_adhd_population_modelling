library(tidyverse)
library(R.matlab)
library(rstatix)
library(perm)
library("dplyr")
library(effsize)
source('util.R')
#Specify parameters of permutation
perm.opt <- permControl(nmc=5000,seed=123451, setSEED = TRUE)

#For the optimal clustering -> calculates significance with respect to control!
thresh <- 0.01

  
#Load .csv used in clustering
data_used_dir <- file.path('data/ASD_ADHD_global.csv')# '/your/path/here/' # <-- EDIT THIS PATH
data_used <- read_csv(data_used_dir)

#Specify file where HYDRA output (in format .mat) is found
results_dir <- file.path('results/sept_results/ASD_ADHD_global_sept/') #"/your/path/here/"  # <-- EDIT THIS PATH
method = 'umap'

#Directory to save results
main_dir <- file.path(results_dir, paste0('FDRcorrected_perm2025_', method))
if (!file.exists(main_dir)){
  dir.create(main_dir, recursive = TRUE)
}


######### MAIN SCRIPT

#Combine all clinical features with data
combined_data <- append_all_features(data_used)

#Select optimal cluster to add to combined data - optimal cluster is saved in csv as cluster_umap or cluster_tsne
optimal_cluster = paste0('cluster_', method)

#Open mat file in results directory
loaded_clustered_data<- read_csv(file.path(results_dir, paste0(method,'_medoids'), 'run_all', 'control_FALSE', '1', paste0(method, '_medoids_clustering.csv'))) 

#Merge clustered data to combined data via the ID column
combined_data <- merge(combined_data, loaded_clustered_data, by = 'ID', all.x = TRUE)

#Change column name to cluster_assignment (used later in the code)
names(combined_data) <- sub(optimal_cluster, 'cluster_assignment', names(combined_data))

#List of clusters  
clusters <- levels(factor(combined_data$cluster_assignment))

#Add controls
combined_data$cluster_assignment[is.na(combined_data$cluster_assignment)] <- -1


#Obtain list of features - EDIT THIS depending on what features you want to analyse
feature_columns <- combined_data %>%
  select(contains('.q.wre'))  %>%
  colnames()


#Remove control group from list of clusters
clusters <- clusters[clusters!='-1']  

#Create empty dataframe to store Cohen's d
cohend_df <- data.frame(matrix( nrow = length(feature_columns), ncol = 7))
colnames(cohend_df)<-c('features', 'control_1', 'p_control_1', 'control_2', 'p_control_2', 'cluster1_cluster2', 'p_cluster1_cluster2')
cohend_df[1] <- feature_columns

########START LOOP ###########


#Iterate over clusters and control
for (cluster in clusters){


  #Create directory to save
  cluster_dir <- file.path(main_dir, paste0('control_', cluster))
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
    if (!is.na(stat_test$p.value) & stat_test$p.value < thresh){
      saving_dir <- file.path(cluster_dir, "sig")
    } else{
      saving_dir <- cluster_dir
    }

    #Attach feature name to saving directory
    filename <- file.path(saving_dir, paste(i, ".png", sep =""))

    #Create plot
    ggsave(filename, create.dir = TRUE)

    #Find Cohen's d - The first argument is the reference group, and the second argument is the comparison group

    if (i !='sex'&& i !='site' && i !='dx.original'){

    find_d <- effsize::cohen.d(data_to_analyse[[i]][data_to_analyse$cluster_assignment==cluster], data_to_analyse[[i]][data_to_analyse$cluster_assignment==-1])

    #Find corresponding row and column where d value should be saved
    cohend_df[which(cohend_df$features == i), current_column] <- find_d$estimate
      p_column <- paste0('p_control_', cluster)
      cohend_df[which(cohend_df$features == i), p_column] <- stat_test$p.value
    }

  }


}
if (length(clusters) == 2){

  #Remove controls
  no_control <- combined_data %>%
    filter(dx.original != 'CN')


  #Create directory to save
  cluster_dir <- file.path(main_dir, 'cluster1_cluster2')
  if (!file.exists(cluster_dir)){
    dir.create(cluster_dir)
  }

  current_column <-paste0('cluster1_cluster2')
  for (i in feature_columns){

    #If entire column is NA - skip
    if (all(is.na(no_control[[i]]))){
      next
    }

    #Remove NA in given feature
    data_to_analyse <- no_control[!is.na(no_control[[i]]), ]

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
      find_d <- effsize::cohen.d(data_to_analyse[[i]][data_to_analyse$cluster_assignment==2], data_to_analyse[[i]][data_to_analyse$cluster_assignment==1])
      cohend_df[which(cohend_df$features == i), current_column] <- find_d$estimate

      #Add p value
      cohend_df[which(cohend_df$features == i), 'p_cluster1_cluster2'] <- stat_test$p.value

    #If p value is significant, include 'significant' directory
    if (!is.na(stat_test$p.value) & stat_test$p.value < thresh){
      saving_dir <- file.path(cluster_dir, "sig")
    } else{
      saving_dir <- cluster_dir
    }

    }
    # #Attach feature name to saving directory
    filename <- file.path(saving_dir, paste(i, ".png", sep =""))
    
    #Create plot
    ggsave(filename, create.dir = TRUE)
  }
}

#FDR correction

#Select columns that include p value
p_columns <- grep("^p_", names(cohend_df), value = TRUE)

# Perform FDR correction for each p_column
for (col in p_columns) {
  p_values <- cohend_df[[col]]

  # Perform FDR correction using p.adjust function
  fdr_corrected <- p.adjust(p_values, method = "fdr")

  # Replace the original p-values with FDR corrected ones
  cohend_df[[paste0(col, '_FDR')]] <- fdr_corrected
}


#Save Cohen's d csv
write_csv(cohend_df, file.path(main_dir, 'Cohens_d_test.csv'))




#CLINICAL ANALYSIS - creates table with differences in clinical scores across groups
features_for_table <- c('age', 'sex', 'IQ', 'RBSR6_TOTAL', 'ADOS_CSS', 'SRS_TOTAL_RAW', 'SWAN_ADHD_I_SUB', 'SWAN_ADHD_HI_SUB', 'FSQC', 'site')
if (nlevels(factor(combined_data$dx.original)) > 2) {
  features_for_table <- c(features_for_table, 'dx.original')
}

columns <- 5 + length(clusters)
table_df <- data.frame(matrix( nrow = length(features_for_table), ncol = columns))
cluster_columns <- paste("cluster_", clusters, sep = "")

colnames(table_df)<-c('features', 'control_mean', 'patient_mean', 'control_pvalue', cluster_columns, 'subtypes_pvalue')
table_df[1] <- features_for_table

#Iterate over features to find significance
for (i in features_for_table){
  print(i)
  
  #Remove NA from data
  data_to_use <-  combined_data[complete.cases(combined_data[[i]]), ]
  no_control <- data_to_use%>%
    filter(cluster_assignment !='-1')
  
  data_to_use <-  data_to_use%>%
    mutate(control_v_all = case_when(dx.original != 'CN' ~ 1,
                                     dx.original == 'CN' ~ -1))
  
  #If feature is sex -> save count and not mean
  if (i =='sex' || (i =='dx.original' && nlevels(factor(combined_data$dx.original))>2) ){
    
    control_sex <- table(combined_data[[i]][combined_data$cluster_assignment =='-1'])
    patient_sex <- table(combined_data[[i]][combined_data$cluster_assignment !='-1'])
    table_df[which(table_df$features == i), 'patient_mean'] <- paste(names(patient_sex), patient_sex, collapse = " ")
    table_df[which(table_df$features == i), 'control_mean'] <-  paste(names(control_sex), control_sex, collapse = " ")
    
    stat_test_control <- chisq.test(data_to_use[[i]], data_to_use$control_v_all)
    stat_test_group <-   chisq.test(no_control[[i]], no_control$cluster_assignment)
    
  }else if (i !='sex'&& i !='dx.original' && i!='site'){
    
    
    #Find mean of feature for control group
    control_mean <- mean(data_to_use[[i]][data_to_use$cluster_assignment==-1], na.rm =TRUE)
    control_sd <- sd(data_to_use[[i]][data_to_use$cluster_assignment==-1], na.rm = TRUE)
    
    #Find mean of feature for patient group
    patient_mean <- mean(data_to_use[[i]][data_to_use$cluster_assignment!=-1], na.rm =TRUE)
    patient_sd <- sd(data_to_use[[i]][data_to_use$cluster_assignment!=-1], na.rm = TRUE)
    
    table_df[which(table_df$features == i), 'control_mean'] <- sprintf("%.2f +/- %.2f", control_mean, control_sd)
    table_df[which(table_df$features == i), 'patient_mean'] <- sprintf("%.2f +/- %.2f", patient_mean, patient_sd)
    
    #Find control vs patients
    stat_test_control <-  kruskal.test(data_to_use[[i]] ~ control_v_all, data = data_to_use)
    stat_test_group <- kruskal.test(no_control[[i]] ~ cluster_assignment, data = no_control)
    
  }
  
  
  #Save p values
  table_df[which(table_df$features == i), 'control_pvalue'] <- stat_test_control$p.value #format(signif(stat_test_control$p.value, 4), scientific = TRUE)
  table_df[which(table_df$features == i), 'subtypes_pvalue'] <- stat_test_group$p.value #format(signif(stat_test_group$p.value, 4), scientific = TRUE)
  
  #Fill in cluster means and counts
  for (cluster in clusters){
    
    #If feature is categorical - count
    if (i =='sex' | (i =='dx.original' && nlevels(factor(combined_data$dx.original))>1 )){
      
      cluster_count <- table(combined_data[[i]][combined_data$cluster_assignment ==cluster])
      table_df[which(table_df$features == i), paste0('cluster_', cluster)] <- paste(names(cluster_count), cluster_count, collapse = " ")
      
      
      #If feature is continuous - mean
    }else if (i !='sex'&& i !='dx.original' && i !='site') {
      cluster_mean <- mean(data_to_use[[i]][data_to_use$cluster_assignment==cluster], na.rm =TRUE)
      cluster_sd <- sd(data_to_use[[i]][data_to_use$cluster_assignment==cluster], na.rm = TRUE)
      table_df[which(table_df$features == i), paste0('cluster_', cluster)] <- sprintf("%.2f +/- %.2f", cluster_mean, cluster_sd)
      
    }
  }
  
}


#FDR correction
#Select columns that include p value
p_columns <- grep("_pvalue$", names(table_df), value = TRUE)

# Perform FDR correction for each p_column
for (col in p_columns) {
  p_values <- table_df[[col]]
  
  # Perform FDR correction using p.adjust function
  fdr_corrected <- p.adjust(p_values, method = "fdr")
  
  # Replace the original p-values with FDR corrected ones
  table_df[[paste0(col, '_FDR')]] <- fdr_corrected
}


#Label significance of features depending on threshold 
table_df <- table_df%>%
  mutate(significance_control = case_when (control_pvalue_FDR < thresh ~'significant',
                                           control_pvalue_FDR > thresh~ 'not significant'))%>%
  mutate(significance_subtype = case_when (subtypes_pvalue_FDR < thresh ~'significant',
                                           subtypes_pvalue_FDR > thresh ~ 'not significant'))


write_csv(table_df, file.path(main_dir, 'table_KW.csv'))

