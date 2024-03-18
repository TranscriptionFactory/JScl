library(tidyverse)
library(effsize)

source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/make_bubbleplots.R')

# # all samples + metadata
rna_data = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/x_data.RDS')
rna_metadata = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/metadata.RDS')

# ##############################################
# # special case: low/high losai for LS LoSAI
rna_metadata$losai = rna_metadata$LoSAI.mLoSSI
#
# # subset to only LS patients
rna_metadata_LS = rna_metadata %>% filter(health == "LS")
rna_data_LS = rna_data[rna_metadata_LS$sample_id, ]
#
# # encode losai (just try median)
rna_metadata_LS$normalized_losai = as.numeric(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/y.csv', row.names = 1)[,1])
#
rna_metadata_LS$losai = as.numeric(rna_metadata_LS$losai)
rna_metadata_LS$normalized_losai = ifelse(rna_metadata_LS$normalized_losai > median(rna_metadata_LS$normalized_losai), "High", "Low")
#
ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
  group_names = c("Peds", "Adult"), LF_to_plot = c(1, 4))
#
rna_metadata_LS2 = rna_metadata_LS
rna_metadata_LS2$onset = rna_metadata_LS2$normalized_losai
#
    lfs = get_LF_data(ls_losai$path, ls_losai$LF_to_plot)
    p = make_bubbleplot(rna_data_LS, lfs, rna_metadata_LS2[, c("onset", "health")], ls_losai$path)
#
#
# ##############################################
# # LS Losai with univariate cor
#
# source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/make_univariate_corr_bubbleplot.R')
#
#
# # add y values to metadata
# rna_metadata_LS$normalized_losai = as.numeric(read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI/y.csv', row.names = 1)[,1])
#     # rna_metadata_LS$LoSAI.mLoSSI = as.numeric(rna_metadata_LS$LoSAI.mLoSSI)
# p = make_univariate_corr_bubbleplot(rna_data_LS, lfs, rna_metadata_LS[, c("onset", "health", "normalized_losai")], ls_losai$path)
# ##############################################
# # LS Onset
# ls_onset = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
#   label = "Blue = Increased in Peds, Red = Increased in Adults",
#   group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))
#
#   lfs = get_LF_data(ls_onset$path, ls_onset$LF_to_plot)
#   p = make_bubbleplot(rna_data_LS, lfs, rna_metadata_LS[, c("onset", "health")], ls_onset$path)
#


# for adult LS vs healthy
# want LFs 4/5
# get LF 4/5 genes
adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), LF_to_plot = c(4, 5) )
    lfs = get_LF_data(adult_ls_vs_healthy$path, adult_ls_vs_healthy$LF_to_plot)
    p = make_bubbleplot_all_lsvsh(rna_data, lfs, rna_metadata[, c("health", "onset")], adult_ls_vs_healthy$path)
#
# all_ls_vs_healthy = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
#   label = 'Blue = Increased in Healthy \n Red = Increased in LS',
#   group_names = c("Healthy", "LS"), LF_to_plot = c(1, 4))
#     lfs = get_LF_data(all_ls_vs_healthy$path, all_ls_vs_healthy$LF_to_plot)
#     p = make_bubbleplot_all_lsvsh(rna_data, lfs, rna_metadata[, c("health", "onset")], all_ls_vs_healthy$path)
#
# # #
# # #
# ls_losai = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
#   label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
#   group_names = NULL, LF_to_plot = c(1, 4))
# # #
# ls_onset = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
#   label = "Blue = Increased in Peds, Red = Increased in Adults",
#   group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))
#   lfs = get_LF_data(ls_onset$path, ls_onset$LF_to_plot)
#
#   p = make_bubbleplot_lsonset(rna_data_LS, lfs, rna_metadata_LS[, c("health", "onset")], ls_onset$path)
#
# #
# #
#   # for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset, ls_losai )) {
#   #   lfs = get_LF_data(m$path, m$LF_to_plot)
#   #   rna_metadata = rna_metadata[rna_metadata$onset == "Adult", ]
#   #   rna_data = rna_data[rna_metadata$sample_id, ]
  #   p = make_bubbleplot(rna_data, lfs, rna_metadata[, c("health", "onset")], m$path)
  # }

# if (F) {
#   # similar set up if want to sort by effect size in original data
#   # # order by effect size
#   # df = get_x_and_y(m$path)
#   #
#   # df_x = df[, -1]
#   # df_x = df_x[, lfs]
#   # df_y = df[, 1]
#
#   # cliffs = filter_by_cliffs_deta(df_x, df_y)
#   #
#   # # take the top 10 by effect size
#   #
#   # eff_size_genes = cliffs$name[1:10]
#   #
#   # # now subset big df
# }
