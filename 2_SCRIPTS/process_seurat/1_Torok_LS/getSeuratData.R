library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)



## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/trh107.h5seurat')
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
  dfBig_H <- dfBig_H[,-which(grepl("HSP|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|^MTND|^ATP", colnames(dfBig_H)))]

  dfBig_H<-as.data.frame(t(dfBig_H))
  dfBig_H$rowVar<- RowVar(dfBig_H)
  dfBig_H$rowMedian<-rowMedian(dfBig_H[,-ncol(dfBig_H)])

  dfBig_H<-dfBig_H[dfBig_H$rowMedian>0,]
  dfBig_H<-dfBig_H[order(dfBig_H$rowVar,decreasing=T),]

  varOrder<-head(dfBig_H,30)
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

write.csv(Var.df2, "INPUT_DATA/VarAll.csv")

check_y<-as.data.frame(data@meta.data)

check_y2<-check_y %>%
  group_by(library_id, health) %>%
  summarise()

check_y2$health <-ifelse(check_y2$health == "LS", 0,1)

write.csv(check_y2, "INPUT_DATA/check_yAll.csv", row.names = F)

