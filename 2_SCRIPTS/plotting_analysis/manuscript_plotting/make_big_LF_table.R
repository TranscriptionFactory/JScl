library(tidyverse)


main_results_folder = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS'

results_folders = list.dirs(main_results_folder, recursive = F)

main_models = results_folders[c(2, 4, 5, 8)]

plot_sig_genes_files = list.files(main_models, pattern = "^plotSigGenes_data.RDS", recursive = T, full.names = T)

model_labels = c("Adult LS vs Healthy", "Peds + Adult LS vs Healthy", "LS mLoSSI", "Peds vs Adult LS")
label_recoding = list(
  c("Higher in Healthy", "Higher in Adult LS"),
  c("Higher in Healthy", "Higher in LS"),
  c("Low mLoSSI", "High mLoSSI"),
  c("Higher in Peds LS", "Higher in Adult LS")
)

all_results = data.frame()
for (t in 1:length(plot_sig_genes_files)) {
  res = readRDS(plot_sig_genes_files[t])

  res$Model = model_labels[t]

  res$color = ifelse(res$color == "Blue",
                     label_recoding[[t]][1], label_recoding[[t]][2])

  all_results = rbind.data.frame(all_results, res)

}

all_results = all_results %>% select(-heights)

names(all_results)[c(1, 5, 6)] = c("Gene", "Association", "LF_number")
saveRDS(all_results, '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/LF_table.RDS')

write_tsv(all_results,
          file = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/LF_table.tsv')



