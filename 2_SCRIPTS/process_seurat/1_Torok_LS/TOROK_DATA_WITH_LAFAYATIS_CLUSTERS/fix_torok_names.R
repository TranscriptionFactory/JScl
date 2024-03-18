library(tidyverse)


lafayatis= as.matrix(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/2_DATA/Var50_mtrp.csv',row.names=1))

LS_transfered_x <- readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/TOROK_DATA_WITH_LAFAYATIS_CLUSTERS/LS_transfered_x2.RDS")
LS_x2 = LS_transfered_x
colnames(LS_x2) = stringr::str_replace(colnames(LS_x2), pattern = "C\\.", "c\\.")

overlapping_in_torok = which(colnames(lafayatis) %in% colnames(LS_x2))
# filter out bad cluster names

lafayatis_clusters = unique(stringr::str_extract(colnames(lafayatis), pattern = "c\\.[0-9]+\\."))

torok_transfered_clusters = which(stringr::str_extract(colnames(LS_x2), pattern = "c\\.[0-9]+\\.") %in% lafayatis_clusters)

LS_x_subset = LS_x2[, torok_transfered_clusters]

write.csv(LS_x_subset, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_MRSS_crosspred/x.csv')


# need to also fix
