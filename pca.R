library(tidyverse)
library(broom)
library(corrr)
library(tidytext)
library(ggplot2)
library(recipes)
source("preprocessing_script.R")

#Parameters
patient = 'both'
crop = NULL
min_age = NULL
max_age = NULL
new_data_filename = NULL #"data/TEST"

file <- "data/combat_centile_ID.csv"
patient_data <- preprocess_data(csv_file = file, patient = patient, min_age = min_age, max_age = max_age, crop = crop, saving_dir = new_data_filename )%>%
  select(-c(1))%>%
  select(!c(session, run, age_days, dx, surfholes, FSQC))

global_feature_data<-patient_data %>%
  select( contains('GMV'), WMVTransformed.q.wre, VentriclesTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre)



pca_recipe <-
  # take all variables
  recipe(~ ., data = patient_data) %>% 
  # specify the ID columns (non-numerical)
  update_role(participant, site, study, group, sex, age, diagnostic, IQ, new_role = 'id') %>% 
  # remove missing values
  step_naomit(all_predictors()) %>% 
  # perform the PCA
  step_pca(all_predictors(), id = "pca") %>% 
  # prepares the recipe by estimating the required parameters
  prep()

patient_pca <- 
  pca_recipe %>% 
  tidy(id = "pca") 

patient_pca

pca_recipe %>% 
  tidy(id = "pca", type = "variance") %>% 
  filter(terms == "percent variance") %>% 
  ggplot(aes(x = component, y = value)) + 
  geom_col() + 
  xlim(c(0, 5)) +
  ylab("% of total variance")

#Transforms data to PC
patient_transform_data <- bake(pca_recipe, new_data = NULL)
ggplot(patient_transform_data, aes(x = PC1, y = PC2, color = sex)) +
  geom_point() +
  labs(title = "PCA",
       x = "Principal Component 1",
       y = "Principal Component 2")

###############
#Global feature
global_recipe <-
  # take all variables
  recipe(~ ., data = global_feature_data) %>% 
  # specify the ID columns (non-numerical)
  # perform the PCA
  step_pca(all_predictors(), id = "pca") %>% 
  # prepares the recipe by estimating the required parameters
  prep()

global_pca <- 
  global_recipe %>% 
  tidy(id = "pca") 

global_pca

global_pca %>% 
  tidy(id = "pca", type = "variance") %>% 
  filter(terms == "percent variance") %>% 
  ggplot(aes(x = component, y = value)) + 
  geom_col() + 
  xlim(c(0, 5)) +
  ylab("% of total variance")

#Transforms data to PC
global_transform_data <- bake(global_recipe, new_data = NULL)
ggplot(global_transform_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(title = "PCA",
       x = "Principal Component 1",
       y = "Principal Component 2")







