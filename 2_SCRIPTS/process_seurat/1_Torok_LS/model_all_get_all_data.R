#!/usr/bin/env Rscript

# Model All LS vs healthy - 50 top var genes

num_var_genes = 5000

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)



## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')
#### ALL SAMPLES

# data <- LS.all

sample_id<-unique(data@meta.data$library_id)

sample_id <- as.character(sample_id)


## GET TOP 30 HIGH VARIANCE GENES PER CLUSTER

# define function for row variance
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

rowMedian<- function(x){
  apply(x, 1, median, na.rm=T)
}



clusters<-c(0:max(as.numeric(data$seurat_clusters)-1))
clusters<-as.factor(clusters)


varList<-c()
dfVar<-NULL

# for each cluster
for (i in 1:length(clusters)
){
  matS<-matrix(nrow=dim(data)[1],ncol=length(sample_id))
  # for each patient
  for (j in 1:length(sample_id)){
    vec<-c()
    #first isolate patient group
    tryCatch(
      {
        tmp<-subset(data ,subset = (data@meta.data[["library_id"]] == sample_id[j])&(seurat_clusters == clusters[i]))
        tmp_df <- as.data.frame(tmp@assays[["RNA"]]@data)
        tmp_df$avg<-rowMedian(tmp_df)
        matS[,j]<-tmp_df$avg
      },
      error= function(cond)
      {
        print(j)
        matS[,j]<-c(rep.int(0,nrow(tmp_df)))
        # print("ERROR")
      }#,
      # finally = print("done")
    )
  }
  dfBig_H<-as.data.frame(matS)

  # repeat gene :/
  dfBig_H<-dfBig_H#[-22664,]

  names<-c()
  clstrs<-substr(clusters[i], 1,2)
  names<-paste(clstrs, rownames(tmp_df),sep=".")
  names<-gsub(" ", "", names)
  names<-gsub("-", "", names)
  names<-gsub("/", "", names)
  names<- paste("c", names, sep=".")
  names<-names#[-22664]
  rownames(dfBig_H)<-names
  rownames(dfBig_H)<-toupper(rownames(dfBig_H))

  dfBig_H<-as.data.frame(t(dfBig_H))
  dfBig_H <- dfBig_H[,-which(grepl("^C.\\d+.HSP|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", stringr::str_to_upper(colnames(dfBig_H))))]

  dfBig_H<-as.data.frame(t(dfBig_H))
  dfBig_H$rowVar<- RowVar(dfBig_H)
  dfBig_H$rowMedian<-rowMedian(dfBig_H[,-ncol(dfBig_H)])

  dfBig_H<-dfBig_H[dfBig_H$rowMedian>0,]
  dfBig_H<-dfBig_H[order(dfBig_H$rowVar,decreasing=T),]

  if (nrow(dfBig_H) < num_var_genes) {
    actual_num_var_genes = nrow(dfBig_H)
  } else {
    actual_num_var_genes = num_var_genes
  }
  varOrder<-head(dfBig_H, actual_num_var_genes)
  if (i==1){
    dfVar<-varOrder
  } else {
    dfVar<-rbind(dfVar, varOrder)
  }
  varList<-c(varList, row.names(varOrder))

}

Var.df<-as.data.frame(t(dfVar[,-length(dfVar)]))


Var.df2<-Var.df
Var.df2[is.na(Var.df2)]<-0
Var.df2<-Var.df2[-nrow(Var.df2),]

rownames(Var.df2)<-sample_id

saveRDS(Var.df2, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/x_data.RDS')

# write.csv(Var.df2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/VarAll_RM.csv")
#
# check_y<-as.data.frame(data@meta.data)
#
# saveRDS(check_y, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/metadata.RDS')
# check_y2<-check_y %>%
#   group_by(library_id, health) %>%
#   summarise()
#
# check_y2$health <-ifelse(check_y2$health == "LS", 1, 0)
#
#
#
# y_ordered = check_y2[match(check_y2$library_id, rownames(Var.df2)), ]
#
# df = cbind.data.frame(y_ordered$health, Var.df2)
# saveRDS(df, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/All_LSvsH.RDS')


#
# # reorder the same as Var.df2
#
# write.csv(check_y2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/check_yAll.csv", row.names = F)

