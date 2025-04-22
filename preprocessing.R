library(tidyverse)
library(filenamer)
source("util.R")

#Parameters
patient = 'both'
crop = NULL
min_age = NULL
max_age =  NULL
new_data_filename = NULL 


file <- "data/combat_centile_ID.csv"
patient_data <- preprocess_data(csv_file = file, patient = patient, min_age = min_age, max_age = max_age, crop = crop, saving_dir = new_data_filename)


#Save all data for UMAP
global_features <-patient_data %>%
    select(ID, contains('GMV'), WMVTransformed.q.wre, VentriclesTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre, group)
write.csv(global_features, "data/all_data.csv", row.names = FALSE)

all_features <-patient_data %>%
    select(ID, contains('q.wre'), group)
 write.csv(all_features, "data/ASD_ADHD_all.csv", row.names=FALSE)


#ADD IDS to sustain data
sustain_asd <- read_csv('data/ADHD_sustain_subtypes_combat_3subtype.csv')[-c(1)]

all_features_data <- read_csv("data/combat_centile_ID.csv")[-c(1)]

sustain_asd$ID <- all_features_data$ID[match(sustain_asd$participant, all_features_data$participant)]

all_features_check <- all_features_data %>%
  filter(ID %in% sustain_asd$ID)
test <- sustain_asd%>%
  filter(site != 'UCSD')%>% 
  filter(site!='UCSDnew')

write.csv(sustain_asd, 'results/ASD_global/sustain_clustering.csv', row.names = FALSE)
