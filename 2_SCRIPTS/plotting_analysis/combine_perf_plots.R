library(tidyverse)


# run1 = list(anno = c("All LS vs H", "Adult LS vs H"),
#   df = readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy_cross_prediction/using_all_ls_vs_healthy/sigLFs_remade_A_Z_cross_pred_performance.RDS"))
#
# run2 = list(anno = c("All LS vs H", "Peds LS vs H"),
#   readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/All_LSvH/sigLFs_cross_pred_performance.RDS"))
#
# run_name_select = c("All LS vs H", "Adult LS vs H", "Peds LS vs H")
run1 = readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy_cross_prediction/using_all_ls_vs_healthy/sigLFs_remade_A_Z_cross_pred_performance.RDS")

run2 = readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/All_LSvH/sigLFs_cross_pred_performance.RDS")

run_name_select = c("All LS vs H", "Adult LS vs H", "Peds LS vs H")

run_list = list(names = c("All LS vs H", "Adult LS vs H", "Peds LS vs H"),
                runs = list(run1$orig_perf, run1$new_perf, run2$new_perf),
                evals = c(run1$orig_eval, run1$new_eval, run2$new_eval))

# very similar to perf_ggplot used in crossprd
combine_perf_plots = function(run_list, title, out_path = NULL, plot_name = "combined_perf") {

  make_df_from_ROCR_performance = function(perf) {

    return(list2DF(list(x = unlist(perf@x.values), y = unlist(perf@y.values))))
  }

  # make into dataframe
  df = data.frame()
  name_title = ""
  for (d in 1:length(run_list$runs)){
    tempdf = make_df_from_ROCR_performance(run_list$runs[[d]])
    tempdf$run = run_list$names[d]

    name_title = paste0(name_title, run_list$names[d], ":", run_list$evals[d], "\n")

    df = rbind(df, tempdf)
  }

  df$run = factor(df$run)

  title = ifelse(!is.null(title), title, "ER Cross Prediction")

  plot_title = paste0(title, "\n",
                      name_title)


  df$x = round(df$x, 1)
  pl = ggpubr::ggline(data = df, x = "x", y = "y", color = "run", plot = "l",
                      size = 3.25, ylab = "True Positive Rate",
                      xlab = "False Positive Rate",
                      palette = "aaas", title = plot_title)


  if(!is.null(out_path)) {
    ggplot2::ggsave(plot = pl, filename = paste0(out_path, "/", plot_name, ".pdf"),
                    height = 6, width = 5.75)
  }

  return(pl)
}


combine_perf_plots(run_list, title = "Generalizing SLIDE LFs", out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures',
                   plot_name = "All_LS_vs_H_perf")
