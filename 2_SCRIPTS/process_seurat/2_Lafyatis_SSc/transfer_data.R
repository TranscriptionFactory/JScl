#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)


# # first run
# if ( F ) {
#   # seurat object is 'harmonskin'
  # load("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/harmonskinV6.RData")

#   torok_data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')
#
#   # lafayatis_to_torok = Seurat::FindTransferAnchors(reference = harmonskin,
#   #                                                  query = torok_data)
#   #
#   # # save(lafayatis_to_torok, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/lafayatis_to_torok.RData')
#   #
#   #
#   # torok_to_lafayatis = Seurat::FindTransferAnchors(reference = torok_data,
#   #                                                  query = harmonskin)
#   # save(torok_to_lafayatis, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis.RData')
#
#
#
#
#   # lafayatis_to_torok = load('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/lafayatis_to_torok.RData')
#   load('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis.RData')
#
#
#   # lafayatis_to_torok_predictions = Seurat::TransferData(anchorset = lafayatis_to_torok,
#   #                                                       refdata = harmonskin@active.ident)
#
#   torok_to_lafayatis_predictions = Seurat::TransferData(anchorset = torok_to_lafayatis,
#                                                         refdata = torok_data@active.ident)
#
#   # save(lafayatis_to_torok, lafayatis_to_torok_predictions,
#   #      file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/lafayatis_to_torok2.RData')
#
#   # save(torok_to_lafayatis, torok_to_lafayatis_predictions,
#   #      file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis2.RData')
#   #
#
#
#
#   saveRDS(torok_to_lafayatis_predictions, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis_predictions.RDS')
# }


############## TRANSFERING LAFAYATIS DATA TO BE IN TOROK ANCHOR SPACE

# called LS.all
# load('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/LS_all_Update_2023.Rdata')
# torok_to_lafayatis = Seurat::FindTransferAnchors(reference = LS.all,
#                                                  query = harmonskin)
# save(torok_to_lafayatis, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_32_serclust_to_lafayatis.RData')
#
# torok_to_lafayatis_predictions = Seurat::TransferData(anchorset = torok_to_lafayatis,
#                                                       refdata = LS.all$seurat_clusters)
# saveRDS(torok_to_lafayatis_predictions, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_32serclust_to_lafayatis_predictions.RDS')




############## TRANSFERING TOROK DATA TO BE IN LAFAYATIS ANCHOR SPACE
load('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/LS_all_Update_2023.Rdata')

load("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/harmonskinV6.RData")

lafayatis_to_torok = Seurat::FindTransferAnchors(reference = harmonskin,
                                                 query = LS.all)

save(lafayatis_to_torok, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_anchors_in_lafayatis.RData')

lafayatis_to_torok_predictions = Seurat::TransferData(anchorset = lafayatis_to_torok,
                                                      refdata = harmonskin$seurat_clusters)

saveRDS(lafayatis_to_torok_predictions, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_anchors_in_lafayatis_predictions.RDS')
