library(tidyverse)
library(ggpubr)

out_dir = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures'
adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), run = "Adult LS vs Healthy")

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), run = "Peds + Adult LS vs Healthy")


ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
  group_names = NULL, run = "Peds + Adult LS LoSAI")

ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "Blue = Increased in Peds, Red = Increased in Adults",
  group_names = c("Peds", "Adult"), run = "Peds + Adult LS Age of Onset")


all_performance_data = list()

# for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset)) {
for (m in list(ls_losai)) {
# for (m in list(all_ls_vs_healthy, ls_onset, adult_ls_vs_healthy)) {


  data = readRDS(paste0(m$path, "/pipeline_step5.rds"))

  data$run = m$run

  all_performance_data = rbind(all_performance_data, data)
}






### Summary showing all performance
performance_plot_data = all_performance_data %>% mutate(method =
                                           recode(method,
                                                  lasso = "Lasso",
                                                  lasso_y = "permuted\nLasso",
                                                  plainER = "True",
                                                  plainER_y = "Permuted")) %>% filter(method %in% c("True", "Permuted"))




for(r in unique(performance_plot_data$run)) {

  temp_df = performance_plot_data %>% filter(run == r)

  # pl1 = ggpubr::ggboxplot(temp_df, x = "method", y = "auc",
  #                         # fill = "method", palette = "aaas", scales = "free_x",
  #                         fill = "method", palette = ggpubr::get_palette("aaas", 10)[8:9], scales = "free_x",
  #                         order = c("True", "Permuted", "Lasso", "permuted\nLasso"),
  #                         xlab = "SLIDE",
  #                         width = 0.6,
  #                         font.xtickslab = c(10),
  #                         ncol = 4) + stat_compare_means(label = "p.signif",
  #                                                        comparisons = list(
  #                                                          c("True", "Permuted")),  #, c("Lasso", "permuted\nLasso")),
  #                                                        step.increase = 0.05, vjust = 0.5)

pl1 = ggpubr::ggboxplot(temp_df, x = "method", y = "corr",
                        # fill = "method", palette = "aaas", scales = "free_x",
                        fill = "method", palette = ggpubr::get_palette("aaas", 10)[8:9], scales = "free_x",
                        xlab = "SLIDE",
                        order = c("True", "Permuted", "Lasso", "permuted\nLasso"),
                        width = 0.6,
                        font.xtickslab = c(10),
                        ncol = 4) + stat_compare_means(label = "p.signif",
                                                       comparisons = list(
                                                         c("True", "Permuted")),  #, c("Lasso", "permuted\nLasso")),
                                                       step.increase = 0.05, vjust = 0.5)

  ggsave(plot = pl1, filename = paste0(out_dir, "/", r, "_performance_with_geno.pdf"), device = "pdf",
         width = 5, height = 4.5)


}




# pl1 = ggpubr::ggboxplot(performance_plot_data, x = "method", y = "auc",
#                         fill = "method", palette = "aaas", scales = "free_x",
#                         order = c("ER", "permuted\nER", "Lasso", "permuted\nLasso"),
#                         facet.by = "run", width = 0.6,
#                         font.xtickslab = c(10),
#                         ncol = 4) + stat_compare_means(label = "p.signif",
#                                                        comparisons = list(
#                                                          c("ER", "permuted\nER"), c("Lasso", "permuted\nLasso")),
#                                                        step.increase = 0.05, vjust = 0.5)
#
# ggsave(plot = pl1, filename = paste0(out_dir, "/all_performance_with_geno.pdf"), device = "pdf",
#        width = 10, height = 5)
