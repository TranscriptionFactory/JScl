# compare the significant latent factors in our two LoSAI runs

library(tidyverse)

orig_dir = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI'
new_dir = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1'

orig_er = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/final_delta_0.1_lambda_1.rds')
new_er = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/orig_er_res.RDS')

orig_z = as.matrix(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/z_mat.csv', row.names = 1))
new_z = as.matrix(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/z_mat.csv', row.names = 1))
orig_slide = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/slide_res.RDS')
new_slide = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/SLIDE_res.RDS')

y = as.matrix(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/y.csv', row.names = 1))

orig_sig_A = orig_er$A[, stringr::str_to_upper(orig_slide$marginal_vars)]
new_sig_A = new_er$A[, new_slide$SLIDE_res$marginal_vars]

matched = intersect(rownames(orig_sig_A), rownames(new_sig_A))

orig_sig_A = orig_sig_A[matched, ]
new_sig_A = new_sig_A[matched, ]
orig_LF_nums = c(1,4)
new_LF_nums = c(1,3)

orig_y_pred = orig_z[, stringr::str_to_upper(orig_slide$marginal_vars)] %*%
  orig_er$beta[as.numeric(stringr::str_remove_all(orig_slide$marginal_vars, pattern = "z")), ]

new_y_pred = new_z[, new_slide$SLIDE_res$marginal_vars] %*%
  new_er$beta[new_slide$SLIDE_res$marginal_vars, ]


orig_y_corr = cor(y, orig_y_pred)
new_y_corr = cor(y, new_y_pred)

cosine_similarity <- function(matrix1, matrix2) {

  # Normalize each row in both matrices
  matrix1_normalized <- matrix1 / sqrt(sum(matrix1^2))#sapply(matrix1, function(x) x / sqrt(sum(x^2)))
  matrix2_normalized <- matrix2 / sqrt(sum(matrix2^2))#sapply(matrix2, function(x) x / sqrt(sum(x^2)))

  # Compute the cosine similarity
  similarity <- matrix1_normalized %*% (matrix2_normalized)

  return(similarity)
}

cosine_similarity_mat <- function(matrix1, matrix2) {
  # Ensure matrices have the same number of columns
  if (ncol(matrix1) != ncol(matrix2)) {
    stop("Both matrices must have the same number of columns")
  }

  # Normalize each row in both matrices
  matrix1_normalized <- apply(matrix1, 2, function(x) x / sqrt(sum(x^2)))
  matrix2_normalized <- apply(matrix2, 2, function(x) x / sqrt(sum(x^2)))

  # Compute the cosine similarity
  similarity <- t(matrix1_normalized) %*% (matrix2_normalized)

  return(similarity)
}

cos_mat = cosine_similarity_mat(orig_sig_A, new_sig_A[, c(1:4, 4)])
pheatmap::pheatmap(cos_mat)

# cell_ext_cos = cosine_similarity(orig_sig_A[, orig_LF_nums[1]], new_sig_A[, new_LF_nums[1]])
cell_int_cos = cosine_similarity(orig_sig_A[, orig_LF_nums[2]], new_sig_A[, new_LF_nums[2]])
