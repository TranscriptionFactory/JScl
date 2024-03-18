#!/usr/bin/env Rscript

# run SLIDE CV on the LoSAI LFs that were used to predict MRSS
library(tidyverse)
library(SLIDE)
# source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/recluster_cross_pred.R')

# source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/remakeA_cross_pred.R')

# cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_to_losai/',
#                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_MRSS_crosspred',
#                      load_slide = load_slide,
#                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_MRSS_crosspred',
#                      orig_new_name = c("MRSS on SSc", "MRSS on LS"))


run_SLIDECV = function(yaml_path) {

  er_input = yaml::yaml.load_file(yaml_path)

  for(j in 1:er_input$nreps){

    withCallingHandlers(SLIDE::paperBench(yaml_path, replicate=j))

  }

  pathLists <- list.files(er_input$out_path,recursive = T,pattern = "results")
  perfList <- lapply(paste0(er_input$out_path,pathLists), readRDS)
  perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


  if (er_input$eval_type == "corr") {

    lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "corr", palette = "aaas",
                                       fill = "method" ) +
      ggpubr::stat_compare_means(label = "p.signif")

  } else {

    lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "auc", palette = "aaas",
                                       fill = "method" ) +
      ggpubr::stat_compare_means(label = "p.signif")
  }

  ggplot2::ggsave(plot = lambda_boxplot, filename = paste0(er_input$out_path, "cv_boxplot.png"), height = 6, width = 6)

  saveRDS(perRes,file=paste0(er_input$out_path,"boxplot_data.rds"))

}


losai_crosspred_SLIDE_location = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/recluster_all_LS_LoSAI/delta0.1_lambda1'

yaml_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_for_MRSS_SLIDECV/params.yaml'

run_SLIDECV(yaml_path)


