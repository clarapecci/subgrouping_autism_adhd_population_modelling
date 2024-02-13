library(tidyverse)
library(filenamer)
source("util.R")

#Parameters
patient = 'both'
crop = NULL
min_age = 10
max_age = 22
new_data_filename = NULL 


file <- "data/combat_centile_ID.csv"
patient_data <- preprocess_data(csv_file = file, patient = patient, min_age = min_age, max_age = max_age, crop = crop, saving_dir = new_data_filename)

#Save all data for UMAP
global_features <-patient_data %>%
    select(ID, contains('GMV'), WMVTransformed.q.wre, VentriclesTransformed.q.wre, totalSA2Transformed.q.wre, meanCT2Transformed.q.wre, group)
write.csv(global_features, "data/ADHD_global_features.csv", row.names = FALSE)

# all_features <-patient_data %>%
#    select(ID, contains('q.wre'), group)
# write.csv(all_features, "data/ADHD_all_features.csv", row.names=FALSE)



# #Save covariate data for HYDRA - transform sex to binary variable
# covariate_data <- patient_data %>%
#   select(ID, sex, age, age_days)%>%
#   mutate(sex = case_when(sex =='Female'~ 1,
#                          sex =='Male'~ 0))
#write.csv(all_data, file.path(new_data_filename+'covariate.csv'))

#Only CT data 
# CT_asd_data <- ASD_data %>%
#   select(ID, contains('CT'), group)
# write.csv(CT_asd_data, 'data/ASD_CT data.csv')
