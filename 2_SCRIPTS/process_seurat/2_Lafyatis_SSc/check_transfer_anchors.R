#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)

# seurat object is 'harmonskin'
load("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/harmonskinV6.RData")

data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')



lafayatis_to_torok = Seurat::FindTransferAnchors(reference = harmonskin,
                                                 query = data)

save(lafayatis_to_torok, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/lafayatis_to_torok.RData')



torok_to_lafayatis = Seurat::FindTransferAnchors(reference = data,
                                                 query = harmonskin)
save(torok_to_lafayatis, file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/torok_to_lafayatis.RData')
