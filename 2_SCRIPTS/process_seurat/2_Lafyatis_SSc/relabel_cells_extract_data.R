#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)

# seurat object is 'harmonskin'
load("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/harmonskinV6.RData")

torok_data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')

torok_to_lafayatis_predictions = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis_predictions.RDS')

harmonskin_pred <- AddMetaData(object = harmonskin, metadata = torok_to_lafayatis_predictions)




######## TOROK data in lafayatis - relabl clusters and save
torok_anchors_in_lafayatis_predictions <- readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_anchors_in_lafayatis_predictions.RDS")
ls.all_pred = AddMetaData(object = torok_data, metadata = torok_anchors_in_lafayatis_predictions)

SaveH5Seurat(ls.all_pred, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/TOROK_CONVERTED_TO_LAFAYATIS_CLUST/LS.all_pred.H5Seurat')
