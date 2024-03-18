library(tidyverse)
library(effsize)

source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/gene_boxplot_effect_size.R')

# all samples + metadata
rna_data = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/x_data.RDS')
rna_metadata = readRDS('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/1_DATA/model_all_get_all_data/metadata.RDS')


# for adult LS vs healthy
# want LFs 4/5
# get LF 4/5 genes
adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), LF_to_plot = c(4, 5) )

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), LF_to_plot = c(1, 4))
#
#
ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
  group_names = NULL, LF_to_plot = c(1, 4))
#
ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "Blue = Increased in Peds, Red = Increased in Adults",
  group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))


for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset, ls_losai )) {
  lfs = get_LF_data(m$path, m$LF_to_plot)
  p = make_gene_boxplot_effect_size(rna_data, lfs, rna_metadata[, c("health", "onset")], m$path)
}

if (F) {
  # similar set up if want to sort by effect size in original data
  # # order by effect size
  # df = get_x_and_y(m$path)
  #
  # df_x = df[, -1]
  # df_x = df_x[, lfs]
  # df_y = df[, 1]

  # cliffs = filter_by_cliffs_deta(df_x, df_y)
  #
  # # take the top 10 by effect size
  #
  # eff_size_genes = cliffs$name[1:10]
  #
  # # now subset big df
}
