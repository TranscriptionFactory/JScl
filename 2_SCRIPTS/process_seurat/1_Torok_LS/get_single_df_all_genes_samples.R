# get all genes/samples

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)



## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')
#### ALL SAMPLES

# data <- LS.all
# onset, health, library_id

data <- subset(x = data, subset = onset != "Peds")

sample_id<-unique(data@meta.data$library_id)

sample_id <- as.character(sample_id)

rna_mat = data@assays[["RNA"]]@data


metadata_annotations = data.frame()

metadata_annotations$sample_id = data$meta.data$library_id
metadata_annotations$health = data$meta.data$health

metadata_annotations$onset = data$meta.data$onset


## GET TOP 20 HIGH VARIANCE GENES PER CLUSTER


clusters<-c(0:max(as.numeric(data$seurat_clusters)-1))
clusters<-as.factor(clusters)


rownames(Var.df2)<-sample_id
saveRDS(Var.df2, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/new_clusters_50/Adult_LSvsH_pre.RDS')


# write.csv(Var.df2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/VarAdult_RM.csv")



check_y<-as.data.frame(data@meta.data)

check_y2<-check_y %>%
  group_by(library_id, health) %>%
  summarise()

check_y2$health <-ifelse(check_y2$health == "LS", 1, 0)

y_ordered = check_y2[match(check_y2$library_id, rownames(Var.df2)), ]

df = cbind.data.frame(y_ordered$health, Var.df2)
saveRDS(df, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/Giffin_new_annotations/Adult_LSvsH.RDS')

#
# write.csv(check_y2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/check_yAdult.csv", row.names = F)

