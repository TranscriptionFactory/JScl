library(tidyverse)
library(pheatmap)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


orig_er_result = readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/sig_genes.RDS")
new_er_result = readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/sig_genes.RDS")

df = data.frame()

for(i in 1:length(orig_er_result)) {

  get_col = function(data) {
    num_unique = apply(data, 2, function(x) length(unique(x)))
    unique_ind = data[, which(num_unique == nrow(data))]

    is_char = sapply(unique_ind[1, ], function(x) is.character(x))
    return(names(unique_ind)[is_char])
  }
  lhs = orig_er_result[[i]]
  lhs = lhs[, get_col(lhs)]

  for (j in 1:length(new_er_result)) {
    rhs = new_er_result[[j]]
    rhs = rhs[, get_col(rhs)]

    # temp = c(i, j, jaccard(lhs, rhs))
    temp = c(i, j, jaccard(unlist(stringr::str_split_i(lhs, pattern = "\\.", i = 3)), unlist(stringr::str_split_i(rhs, pattern = "\\.", i = 3))))

    df = rbind.data.frame(temp, df)
  }
}

df_mat = matrix(df[, 3], length(orig_er_result), length(new_er_result))

rownames(df_mat) = length(orig_er_result):1
colnames(df_mat) = length(new_er_result):1


pheatmap(df_mat)
