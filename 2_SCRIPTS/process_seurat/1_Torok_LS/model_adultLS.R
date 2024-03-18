# Model Adult LS vs Healthy - 50 top var genes

num_var_genes = 50

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)



## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')
#### ALL SAMPLES

# data <- LS.all
data <- subset(x = data, subset = (onset == "Adult" & health == "LS"))

sample_id<-unique(data@meta.data$library_id)

sample_id <- as.character(sample_id)


## GET TOP 20 HIGH VARIANCE GENES PER CLUSTER


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

  varOrder<-head(dfBig_H,num_var_genes)
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


# write.csv(Var.df2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/VarAdultLS_RM.csv")
saveRDS(Var.df2, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/new_clusters_50/Adult_LoSAI_pre.RDS')


check_y<-as.data.frame(data@meta.data)

check_y2<-check_y %>%
  group_by(library_id, health) %>%
  summarise()

check_y2$health <-ifelse(check_y2$health == "LS", 1, 0)

get_yvals = function(ymat, metapath = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/scRNAseq_LS_Update_Metadata_Jacob.csv'){
  metadata = read.csv(metapath)

  return(ymat %>% left_join(metadata, by = c("library_id" = "sample_id")))
}

check_y2 = get_yvals(check_y2)

check_y2$LoSAI.mLoSSI = as.numeric(check_y2$LoSAI.mLoSSI)

resid = check_y2 %>% select(sex, LoSAI.mLoSSI)

resid$sex = ifelse(resid$sex == "F", 1, 0)
rfit = lm(LoSAI.mLoSSI ~ ., data = resid)

check_y2$LoSAI.mLoSSI = rfit$residuals


check_y2 = check_y2[, c("library_id", "LoSAI.mLoSSI")]


y_ordered = check_y2[match(check_y2$library_id, rownames(Var.df2)), ]

df = cbind.data.frame(y_ordered$LoSAI.mLoSSI, Var.df2)
saveRDS(df, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/Giffin_new_annotations/Adult_LoSAI.RDS')


# write.csv(check_y2, "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/data/check_yAdult.csv", row.names = F)
#
