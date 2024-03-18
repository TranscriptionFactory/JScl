library(tidyverse)
library(caret)
library(ggpubr)


y = as.matrix(read.csv(
  '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/orig_y.csv',
  row.names = 1))

z = as.matrix(read.csv(
  '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/orig_z.csv',
  row.names = 1))

slide_res = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/SLIDE_res.RDS')



z = z[, slide_res$SLIDE_res$marginal_vars]

num_k_folds = 4
num_runs = 10


for (run_num in 1:num_runs) {
  k_folds = caret::createFolds(y, k = num_k_folds)

  for(i in 1:num_k_folds) {

    train_x = z[[]]
  }

