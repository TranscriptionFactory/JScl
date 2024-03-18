library(tidyverse)


rna_data = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/x_data.RDS')
rna_metadata = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/metadata.RDS')



# peds samples
peds_data = rna_data[which(rna_metadata$onset == "Peds"), ]
peds_metadata = rna_metadata[which(rna_metadata$onset == "Peds"), ]
y = peds_metadata$health
y = ifelse(y == "LS", 1, 0)
write.csv(peds_data, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/x.csv')
write.csv(y, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/y.csv')
