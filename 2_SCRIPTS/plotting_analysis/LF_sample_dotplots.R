library(tidyverse)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(qgraph) ## for making the network


get_files_from_internal_run_path = function(run_path,
                                            file_names = list(x = "x.csv",
                                                              y = "y.csv",
                                                              slide_res = "slide_res.RDS",
                                                              # lasso_res = "lasso_results.RDS",
                                                              # sig_genes = "plotSigGenes_data.RDS",
                                                              z = "z_mat.csv")) {

  check_num_file_matches = function(file_path, pattern) {

    file_matches = list.files(file_path, full.names = T, pattern = pattern, recursive = T)

    # file_paths = list.files(file_path, full.names = T, recursive = T)
    #
    # file_matches = file_paths[stringr::str_which(file_paths, pattern = pattern)]

    if (length(file_matches) > 0) {
      return(file_matches[1])
    }
    return()
  }

  # main_run_dir = stringr::str_split_i(run_path, pattern = "/pipeline3", i = 1)

  data = file_names

  for (i in 1:length(file_names)) {

    # file_path = check_num_file_matches(main_run_dir, file_names[[i]])
    file_path = check_num_file_matches(run_path, file_names[[i]])

    if (stringr::str_detect(file_names[[i]], pattern = ".csv")) {

      data[[i]] = read.csv(file_path, row.names = 1)

    } else if (stringr::str_detect(file_names[[i]], pattern = ".RDS")) {
      data[[i]] = readRDS(file_path)
    }
  }
  return(data)
}


make_LF_sample_dotplot = function(df, out_path, file_name, plot_height = 10, plot_width = 15) {
  # df = [y, Z1, Z2]

  pl = ggpubr::ggscatter(df, x = colnames(df)[2], y = colnames(df)[3], color = colnames(df)[1],size = 7) +
    scale_color_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) +
  theme(legend.position = "left")

  ggplot2::ggsave(plot = pl, filename = paste0(out_path, "/", file_name, ".png"),
                  height = 4.5, width = 7)
}

getChullPolygon = function(data) {
  results = data.frame()
  for (group in unique(data[,1])) {
    # filter each group
    df = data[data[,1] == group, ]
    df[,2] = as.numeric(df[,2])
    df[,3] = as.numeric(df[,3])

    # df = data %>% dplyr::filter(colnames(data)[1] == group)
    outliers1 = graphics::boxplot(df[,2], plot = F)$out
    outliers2 = graphics::boxplot(df[,3], plot = F)$out

    chdf = df %>% dplyr::filter(!df[,2] %in% c(outliers1) & !df[,3] %in% c(outliers2))

    boundary = grDevices::chull(chdf)
    BumpX = chdf[,2][boundary] #+ 0.1*(df$Comp1[boundary] - mx)
    BumpY = chdf[,3][boundary] #+ 0.1*(df$Comp2[boundary] - my)

    results = rbind(results, list(Comp1 = BumpX, Comp2 = BumpY, True = chdf[,1][boundary]))
  }
  return(results)
}

make_LF_sample_dotplot_chull = function(df, out_path, file_name, plot_height = 10, plot_width = 15) {
  # df = [y, Z1, Z2]

  chull_df = getChullPolygon(df)


  pl = ggpubr::ggscatter(df, x = colnames(df)[2], y = colnames(df)[3], color = colnames(df)[1],size = 7) +
     # geom_polygon(data = data.frame(chull_df), aes(x = Comp1, y = Comp2, fill = True), alpha = 0.25, inherit.aes = T) +
    stat_chull(aes(fill = as.factor(df[,1])), alpha = 0.2, geom = "polygon", show.legend = F) +
    scale_color_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) +
    scale_fill_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) +

    # scale_fill_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) +
    theme(legend.position = "left")

  ggplot2::ggsave(plot = pl, filename = paste0(out_path, "/chull_", file_name, ".png"),
                  height = 4.5, width = 7)
}


plot_LF_sample_dotplot = function(run_path, encodings, lf_nums, out_folder = "/LF_sample_dotplots/") {

  out_path = paste0(run_path, out_folder)
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = T)
  }
  data = get_files_from_internal_run_path(run_path)

  data$y[,1] = ifelse(data$y[,1] == 0, encodings[1], encodings[2])

  data$y[, 1] = factor(data$y[,1], levels = encodings)

  data$df = cbind.data.frame(data$y, data$x)
  z_cols = stringr::str_to_upper(data$slide_res$marginal_vars)[lf_nums]

  make_LF_sample_dotplot_chull(cbind.data.frame(data$y, data$z[, z_cols]), out_path, paste0("LF", "_", lf_nums[1], "_", lf_nums[2]))
}



adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  group_names = c("Healthy", "LS"), LF_to_plot = c(4, 5) )

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  group_names = c("Healthy", "LS"), LF_to_plot = c(1, 4))

# ls_losai = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
#   group_names = , LF_to_plot = c(1, 4))
# #
ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))



for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset)) {
  plot_LF_sample_dotplot(m$path, m$group_names, m$LF_to_plot)
}



