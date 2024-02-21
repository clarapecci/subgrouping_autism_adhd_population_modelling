library(tidyverse)
library(broom)
library(corrr)
library(tidytext)
library(ggplot2)
library(recipes)
source("util.R")
library(FactoMineR)


#1.1 Generate data from original csv
#Parameters
# patient = 'both'
# crop = NULL
# min_age = NULL
# max_age = NULL
# new_data_filename = NULL #"data/TEST"
# 
# file <- "data/combat_centile_ID.csv"
# patient_data <- preprocess_data(csv_file = file, patient = patient, min_age = min_age, max_age = max_age, crop = crop, saving_dir = new_data_filename )%>%
#   select(-c(1))%>%
#   select(!c(session, run, age_days, dx, surfholes, FSQC))

#1.2 Open csv directly
patient_data <-read_csv("data/ASD_all_features.csv")

#1.3 Saving directory 
saving_dir <- 'figures/ASD_all_pca'

#2. Select global features 
global_feature_data<-patient_data %>%
  select( contains('GMV'), WMVTransformed.q.wre, VentriclesTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre)

###### COMBINE FEATURES IN SINGLE DATAFRAME
#Load all clinical data + original features
clinical_data <-read_csv("data/clinical_data_ID.csv")[-c(1)]
all_features_data <- read_csv("data/combat_centile_ID.csv")[-c(1)]

#Check IDs used in clustering to extract desired features from original + clinical data
used_ID <- patient_data$ID

#Extract correct clinical features
clinical_features <- clinical_data %>% 
  filter(ID %in% used_ID) %>%
  select (-c(participant, site, dx.original, sex, ID))

#Extract non clinical features not used in clustering
non_clinical_features <- all_features_data %>% 
  filter(ID %in% used_ID) %>%
  select(site, age, sex, IQ, dx.original, dx.original)

#Combine clinical and non clinical features
patient_data <- patient_data %>%
  cbind(clinical_features)%>%
  cbind(non_clinical_features)


###############

fill ='site'

pca_result <- PCA(global_feature_data, scale.unit = FALSE)

patient_data_combined <- patient_data%>%
  mutate(PC1 = pca_result$ind$coord[, 1])%>%
  mutate(PC2 = pca_result$ind$coord[, 2])
  
ggplot(patient_data_combined, aes(x = PC1, y = PC2, color = !!sym(fill))) +
  geom_point() + #scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Scatter Plot", x = "PC1", y = "PC2")

ggsave(file.path(paste0(saving_dir, '_', fill, '.png')))


#Find which features carry most variance per principal component
# Extract loadings for the first principal component
loadings_pc1 <- pca_result$var$coord[, 2]

# Identify features with the highest absolute loadings
top_features_pc1 <- names(sort(abs(loadings_pc1), decreasing = TRUE))

# Display the top features for the first principal component
print(top_features_pc1)





