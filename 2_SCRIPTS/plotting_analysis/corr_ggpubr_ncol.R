# plot performance plots in columns instead of rows
library(tidyverse)
library(ggpubr)

cor_ggpubr_col = function(orig_perf, new_perf, title = NULL, orig_new_name = NULL) {

  if (is.null(orig_new_name)) {
    orig_new_name = c("ER", "crosspred")
  }

  # combine dfs
  orig_perf$run = orig_new_name[1]
  new_perf$run = orig_new_name[2]


  df = rbind.data.frame(orig_perf, new_perf)

  title = ifelse(!is.null(title), title, "ER Cross Prediction")

  pl = ggpubr:: ggscatter(df, x = "true",
                          y = "pred",
                          palette = "npg",
                          xlab = "True Y",
                          ylab = "Predicted Y",
                          facet.by = "run",
                          ncol= 2,
                          scales = "free",
                          add = "reg.line",
                          title = title,
                          add.params = list(color = "blue", fill = "lightgray"),
                          conf.int = T) +
    stat_cor(method = "spearman")

  return(pl)
}



perf_data <- readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1/sigLFs_cross_pred_performance.RDS")

output_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1'


perf_data$orig_perf$pred = perf_data$orig_perf$pred/max(perf_data$orig_perf$pred)#scale(perf_data$orig_perf$pred, T, F)
perf_data$new_perf$pred = perf_data$new_perf$pred/max(perf_data$new_perf$pred)#scale(perf_data$new_perf$pred, T, F)
pl = cor_ggpubr_col(perf_data$orig_perf, perf_data$new_perf, title = "",
                    orig_new_name = c("LS mLoSSI", "SSc MRSS"))




ggplot2::ggsave(filename = paste0(output_path, "/col_crosspred_perf_plot.pdf"), plot = pl,
                height = 4, width = 9)



