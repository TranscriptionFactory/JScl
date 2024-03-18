library(tidyverse)
library(ggpubr)

# Boxplots showing latent factor values for each sample. Make sure you have
# X and Y dataframes in the results folder, or nearby.


plot_LF_boxplots = function(path, df_filename_RDS_pattern = "df.RDS",
                            group_names = NULL) {


  check_for_file = function(path, file_pattern) {

    f = list.files(path, recursive = T, full.names = T, pattern = file_pattern)
    if( length(f) > 0 ) {
      return(f[1])
    } else {
      cat("\n Folder must have exactly one file with desired name: ", pattern, "\n")
      return(NULL)
    }
  }

  get_x_and_y = function(path) {
    x = check_for_file(path, file_pattern = "x.csv")
    y = check_for_file(path, file_pattern = "y.csv")

    load_matrix_csv = function(path) {
      return(read.csv(path, row.names = 1))
    }

    if (all(!is.null(x), !is.null(y))) {


      reorder_y_to_xrows = function(x, y) {

        if (all(rownames(x) %in% rownames(y))) {
          reordered_y = as.matrix(y)[rownames(x), ]
          # return(as.matrix(reordered_y))
          return(reordered_y)
        }
        return(y)
      }

      x = load_matrix_csv(x)
      y = load_matrix_csv(y)

      y = reorder_y_to_xrows(x, y)

      return(cbind.data.frame(y, x))


      # return(cbind.data.frame(load_matrix_csv(y), load_matrix_csv(x)))
    } else {
      df = check_for_file(path, file_pattern = df_filename_RDS_pattern)
      if (!is.null(df)) {
        return(readRDS(df))
      } else {
        return(NULL)
      }
    }
  }


  df = get_x_and_y(path)
  if (is.null(df)) {
    cat("\n Failed to load dataframes. Check path \n")
    return()
  } else {
    names(df)[1] = "y"
  }


  # load Z matrix
  z_mat = check_for_file(path, file_pattern = "z_mat.csv")

  if ( is.null(z_mat)) {
    cat("\n Z matrix not found. Check path \n")
    return()
  }

  sig_genes_data = check_for_file(path, file_pattern = "slide_res.RDS")

  all_LFs = NULL

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(unique(readRDS(sig_genes_data)$marginal_vars))

  } else {
    cat("\n Couldn't find significant latent factors. Plotting all latent factors \n")
    all_LFs = colnames(z_mat)
  }


  z_mat = as.matrix(read.csv(z_mat, row.names = 1))[, all_LFs]


  # Have X, Y and Z matrix loaded. Can plot latent factors now

  # first attach Y to Z matrix
  z_combined = cbind.data.frame(df$y, z_mat)

  names(z_combined)[1] = "y"

  z_pivot = z_combined %>% tidyr::pivot_longer(.,
                                               cols = colnames(z_mat),
                                               names_to = "LF_num",
                                               values_to = "LF_val")

  # rename latent factors
  for (z in 1:length(all_LFs)) {

    z_pivot$LF_num = ifelse(z_pivot$LF_num == all_LFs[z], z, z_pivot$LF_num)
  }

  z_pivot$LF_val = as.numeric(z_pivot$LF_val)
  z_pivot$y = as.numeric(z_pivot$y)

  if (!is.null(group_names) && length(group_names) >= 2) {
    z_pivot$y = ifelse(z_pivot$y == min(z_pivot$y), group_names[1], group_names[2])

    z_pivot$y = factor(z_pivot$y)

  }

  dir_name = paste0(path, "/LF_boxplots")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }
  original_wd = getwd()

  for(LF in unique(z_pivot$LF_num)) {
    LF_temp = z_pivot %>% filter(LF_num == LF)


    # plots for regression case
    is_classification = length(unique(z_pivot$y)) == 2

    if (is_classification) {


      # if (!is.null(group_names) && length(group_names) >= 2) {
      #   LF_temp$y = ifelse(LF_temp$y == min(LF_temp$y), group_names[1], group_names[2])
      # }

      pl = ggpubr::ggboxplot(LF_temp,
                             x = "y",
                             y = "LF_val", palette = "aaas",
                             fill = "y",
                             xlab = " ",
                             ylab = paste0("Latent Factor #", LF) ) +
        ggpubr::stat_compare_means(method = "wilcox.test",
                                   symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                      symbols = c("*", "*", "*", "*", ".")),
                                   aes(label = after_stat(p.signif)))

      plot_title = paste0(dir_name, "/LF", LF, "_boxplot.pdf")
      ggplot2::ggsave(filename = plot_title,
                      plot = pl, height = 3.5, width = 4,
                      device = "pdf")

      # save all LF plots in one

      pl = ggpubr::ggboxplot(z_pivot,
                             x = "y",
                             y = "LF_val",
                             facet.by = "LF_num",
                             palette = "aaas",
                             fill = "y",
                             xlab = "Group",
                             ylab = "Latent Factor Value" ) +
        ggpubr::stat_compare_means(method = "wilcox.test")

      plot_title = paste0(dir_name, "/", "all_boxplot.pdf")
      ggplot2::ggsave(filename = plot_title,
                      plot = pl, height = 10, width = 10,
                      device = "pdf")


    } else {

      # sort by y value
      # LF_temp = LF_temp[order(LF_temp$y), ]
      LF_temp$y = round(as.numeric(LF_temp$y), 3)



      pl = ggpubr::ggscatter(LF_temp, x = "y",
                              y = "LF_val",
                              palette = "npg",
                              xlab = "LoSAI",
                              add = "reg.line",
                              add.params = list(color = "blue", fill = "lightgray"),
                             conf.int = T,
                              ylab = paste0("Latent Factor #", LF) ) +
        stat_cor(method = "pearson")

      plot_title = paste0(dir_name, "/LF", LF, "_regplot.pdf")
      ggplot2::ggsave(filename = plot_title,
                      plot = pl, height = 3.5, width = 4,
                      device = "pdf")

      pl = ggpubr:: ggscatter(z_pivot, x = "y",
                              y = "LF_val",
                              palette = "npg",
                              xlab = "LoSAI",
                              facet.by = "LF_num",
                              add = "reg.line",
                              add.params = list(color = "blue", fill = "lightgray"),
                              conf.int = T,
                              ylab = "Significant Latent Factors" ) +
        stat_cor(method = "pearson")

      plot_title = paste0(dir_name, "/", "all_regplot.pdf")
      ggplot2::ggsave(filename = plot_title,
                      plot = pl, height = 10, width = 10,
                      device = "pdf")

    }
  }
  setwd(original_wd)

}


adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"))

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"))


ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
  group_names = NULL)

ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "Blue = Increased in Peds, Red = Increased in Adults",
  group_names = c("Peds", "Adult"))

for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_losai, ls_onset)) {
# for (m in list(ls_losai)) {

  plot_LF_boxplots(m$path, group_names = m$group_names)
}




