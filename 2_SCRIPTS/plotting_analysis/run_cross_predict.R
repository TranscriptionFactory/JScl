source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/remakeA_cross_pred.R')

# load_slide = F


for (load_slide in c(T, F)) {


  ########## predict losai
  cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_to_losai/',
                       new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_MRSS_crosspred',
                       load_slide = load_slide,
                       out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI_MRSS_crosspred',
                       orig_new_name = c("MRSS on SSc", "MRSS on LS"))

  ########## losai to mrss
  cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
                       new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss',
                       load_slide = load_slide,
                       out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_use_losai_to_mrss/using_all_LS_LoSAI',
                       orig_new_name = c("LoSAI on LS", "LoSAI on SSc"))

  ########## LS all to LS Adult
  cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
                       new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy_cross_prediction',
                       load_slide = load_slide,
                       out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy_cross_prediction/using_all_ls_vs_healthy/',
                       orig_new_name = c("All LS vs H", "Adult LS vs H"))




  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/Adult_LSvH',
  #                      orig_new_name = c("Adult LS vs H", "Peds LS vs H"))
  #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/All_LSvH',
  #                      orig_new_name = c("All LS vs H", "Peds LS vs H"))
  #



  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred',
  #                      orig_new_name = c("All LS vs H", "SSc vs H"))
  #
  #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred/using_Adult_LS_vs_H/',
  #                      orig_new_name = c("Adult LS vs H", "SSc vs H"))
  #
  #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_no_tofa',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_no_tofa/using_All_LS_vs_H',
  #                      orig_new_name = c("All LS vs H", "SSc vs H"))
  #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_no_tofa',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/6_LAFAYATIS_DATA/cross_pred_no_tofa/using_Adult_LS_vs_H/',
  #                      orig_new_name = c("Adult LS vs H", "SSc vs H"))
  #
  #
  # # # onset
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/Adult_LSvH',
  #                      orig_new_name = c("Adult LS vs H", "Peds LS vs H"))
  # #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
  #                      load_slide = load_slide,
  #                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/All_LSvH',
  #                      orig_new_name = c("All LS vs H", "Peds LS vs H"))
  # #
  # # # losai
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LoSAI_cross_prediction',
  #                      load_slide = load_slide,
  #                      orig_new_name = c("All LS LoSAI", "Adult LS LoSAI"))
  # #
  # #
  # cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  #                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_LoSAI_cross_prediction',
  #                      load_slide = load_slide,
  #                      orig_new_name = c("All LS LoSAI", "Peds LS LoSAI"))



}




# source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/cross_predict.R')
#
# load_slide = F
#
# cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
#                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
#                      load_slide = load_slide,
#                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/Adult_LSvH',
#                      orig_new_name = c("Adult LS vs H", "Peds LS vs H"))
#
# cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
#                      new_run_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction',
#                      load_slide = load_slide,
#                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/All_LSvH',
#                      orig_new_name = c("All LS vs H", "Peds LS vs H"))
#
#
#



# cross_predict_loader(orig_run_path = '/ix/djishnu/Aaron/3_collabs/2_Eickelberg_SSc/longRunData/ER_NoRM/results/SSCILD_COPD/pipeline3/delta0.4_lambda1',
#                      new_run_path = '/ix/djishnu/Aaron/3_collabs/2_Eickelberg_SSc/longRunData/ER_NoRM/results/IPF_SSCILD',
#                      load_slide = F,
#                      out_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Peds_cross_prediction/eb_test',
#                      orig_new_name = c("SSCILD vs COPD", "IPF vs SSCILD"))
