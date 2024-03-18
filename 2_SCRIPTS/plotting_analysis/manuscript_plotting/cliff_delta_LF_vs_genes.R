library(effsize)
library(tidyverse)
library(ggpubr)

cliff_delta_manual = function(data, y) {
  results = matrix(data = 0, nrow = ncol(data), ncol = 1, dimnames = list(colnames(data), "res"))

  for (col in colnames(data)) {

    cd = cliff.delta(data[y == 0, col], data[y == 1, col], use.unbiased = F)$estimate
    # cd = cliff.delta(d = data[, col], f = y)$estimate

    results[col, "res"] = cd

  }
  return(results)
}



plot_LF_boxplots = function(path, df_filename_RDS_pattern = "df.RDS",
                            group_names = NULL, LF_to_plot = NULL) {

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
      return(as.matrix(read.csv(path, row.names = 1)))
    }

    if (all(!is.null(x), !is.null(y))) {


      reorder_y_to_xrows = function(x, y) {

        y = as.matrix(y[match(rownames(x), rownames(y)),])
        rownames(y) = rownames(x)
        return(y)
      }

      x = load_matrix_csv(x)
      y = load_matrix_csv(y)

      y = reorder_y_to_xrows(x, y)

      return(cbind.data.frame(y, x))

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

  # load slide results
  sig_genes_data = readRDS(check_for_file(path, file_pattern = "sig_genes.RDS"))

  plotted_sig_genes = readRDS(check_for_file(path, file_pattern = "plotSigGenes_data.RDS"))

  # load Z matrix
  z_mat_path = check_for_file(path, file_pattern = "z_mat.csv")

  if ( is.null(z_mat_path)) {
    cat("\n Z matrix not found. Check path \n")
    return()
  }

  z_mat_full = as.matrix(read.csv(z_mat_path, row.names = 1))


  slide_res_path = check_for_file(path, file_pattern = "slide_res.RDS")

  all_LFs = NULL

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(unique(readRDS(slide_res_path)$marginal_vars))

  } else {
    cat("\n Couldn't find significant latent factors. Plotting all latent factors \n")
    all_LFs = colnames(z_mat_full)
  }


  dir_name = paste0(path, "/LF_cliffs_delta")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }
  original_wd = getwd()

  if (!is.null(LF_to_plot)) {
    all_LFs = all_LFs[LF_to_plot]
  }

  z_mat = z_mat_full[, all_LFs]

  z_combined = cbind.data.frame(df$y, z_mat)
  names(z_combined)[1] = "y"


  # get gene data
  sig_genes_vec = unlist(lapply(sig_genes_data, function(x) x['gene']))

  sig_genes_plotted_vec = plotted_sig_genes$name

  # gene_mat = df[, which(!colnames(df) %in% sig_genes_vec)]
  gene_mat = df[, setdiff(colnames(df), sig_genes_vec)]

  sample_size = length(sig_genes_plotted_vec)

  sampled_non_sig_genes = sample(2:ncol(gene_mat), sample_size)
  # sample DF
  # gene_mat = df[, sampled_non_sig_genes]

  gene_mat = gene_mat[, sampled_non_sig_genes]

  # add y back
  gene_mat = cbind.data.frame(df$y, gene_mat)
  names(gene_mat)[1] = "y"


  # run cliffs delta

  doing_classification = ifelse(length(unique(df$y)) == 2, T, F)

  if (doing_classification) {

    z_combined$y = as.numeric(z_combined$y)

    # cliff_slide = abs(apply(z_combined[, -1], 2, function(x) cliff.delta(d = x, f = z_combined$y, return.dm = T)$estimate))
    # cliff_slide = abs(apply(z_combined[, -1], 2, function(x) cliff.delta(x[z_combined$y == 0], x[z_combined$y == 1], return.dm = T)$estimate))
    cliff_slide = abs(cliff_delta_manual(z_combined[, -1], z_combined[, 1]))


    cliff_df = data.frame(cliff_slide)
    colnames(cliff_df)[1] = "cliff"
    cliff_df$source = paste0("SZ", 1:length(rownames(cliff_df)))
    cliff_df$condition = "SLIDE"

    # cliff_genes = abs(apply(gene_mat[, -1], 2, function(x) cliff.delta(d = x, f = gene_mat$y, return.dm = T)$estimate))
    # cliff_genes = abs(apply(gene_mat[, -1], 2, function(x) cliff.delta(x[gene_mat$y == 0], x[gene_mat$y == 1], return.dm = T)$estimate))

    cliff_genes = abs(cliff_delta_manual(gene_mat[, -1], gene_mat[, 1]))
    cliff_df2 = data.frame(cliff_genes)
    colnames(cliff_df2)[1] = "cliff"
    cliff_df2$source = ""#colnames(gene_mat[, -1])
    cliff_df2$condition = "Nonsig Genes"

    plot_lab = "Cliffs Delta"

  } else {
    cliff_slide = apply(z_combined[, -1], 2, function(x) abs(cor(x, z_combined$y)))
    cliff_df = data.frame(cliff_slide)
    colnames(cliff_df)[1] = "cliff"
    cliff_df$source = paste0("SZ", 1:length(rownames(cliff_df)))
    cliff_df$condition = "SLIDE"

    cliff_genes = apply(gene_mat[, -1], 2, function(x) abs(cor(x, gene_mat$y)))
    cliff_df2 = data.frame(cliff_genes)
    colnames(cliff_df2)[1] = "cliff"
    cliff_df2$source = ""#colnames(gene_mat[, -1])
    cliff_df2$condition = "Nonsig Genes"

    plot_lab = "Correlation with Y"

  }


  cliff_combined = rbind(cliff_df, cliff_df2)

  saveRDS(list(cliffs_delta = cliff_combined,
               sampled_non_sig_genes = sampled_non_sig_genes),
          paste0(dir_name, "/cliff_results.RDS"))


  pl = ggpubr::ggboxplot(cliff_combined, x = "condition", y = "cliff",
                         xlab = "Group", repel = T, fill = "condition",
                         palette = "aaas",
                         ylab = plot_lab, add = "point", label = "source") +
    ggpubr::stat_compare_means(method = "wilcox.test", label.x.npc = "right")


  plot_title = paste0(dir_name, "/CD_boxplot_LF_vs_nonsig.png")

  ggplot2::ggsave(filename = plot_title,
                  plot = pl, height = 7, width = 7,
                  device = "png")

}


adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), LF_to_plot = c(4, 5) )

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \n Red = Increased in LS',
  group_names = c("Healthy", "LS"), LF_to_plot = c(1, 4))


ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
  group_names = NULL, LF_to_plot = c(1, 4))

ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "Blue = Increased in Peds, Red = Increased in Adults",
  group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))

for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_losai, ls_onset)) {
  plot_LF_boxplots(m$path, group_names = m$group_names)
}













#
# plot_LF_boxplots = function(path, df_filename_RDS_pattern = "df.RDS",
#                             group_names = NULL, LF_to_plot = NULL) {
#
#
#   check_for_file = function(path, file_pattern) {
#
#     f = list.files(path, recursive = T, full.names = T, pattern = file_pattern)
#     if( length(f) > 0 ) {
#       return(f[1])
#     } else {
#       cat("\n Folder must have exactly one file with desired name: ", pattern, "\n")
#       return(NULL)
#     }
#   }
#
#   get_x_and_y = function(path) {
#     x = check_for_file(path, file_pattern = "x.csv")
#     y = check_for_file(path, file_pattern = "y.csv")
#
#     load_matrix_csv = function(path) {
#       return(as.matrix(read.csv(path, row.names = 1)))
#     }
#
#     if (all(!is.null(x), !is.null(y))) {
#
#
#       reorder_y_to_xrows = function(x, y) {
#
#         # y = as.matrix(y)
#         # if (all(rownames(x) %in% rownames(y))) {
#         #   reordered_y = as.matrix(y)[rownames(x), ]
#         #   # return(as.matrix(reordered_y))
#         # return(reordered_y)
#         return(as.matrix(y)[match(rownames(x), rownames(y))])
#
#         # } else {
#         # return(NULL)
#         # }
#       }
#
#       x = load_matrix_csv(x)
#       y = load_matrix_csv(y)
#
#       y = reorder_y_to_xrows(x, y)
#
#       return(cbind.data.frame(y, x))
#
#
#       # return(cbind.data.frame(load_matrix_csv(y), load_matrix_csv(x)))
#     } else {
#       df = check_for_file(path, file_pattern = df_filename_RDS_pattern)
#       if (!is.null(df)) {
#         return(readRDS(df))
#       } else {
#         return(NULL)
#       }
#     }
#   }
#
#
#   df = get_x_and_y(path)
#   if (is.null(df)) {
#     cat("\n Failed to load dataframes. Check path \n")
#     return()
#   } else {
#     names(df)[1] = "y"
#   }
#
#   # load slide results
#   sig_genes_data = readRDS(check_for_file(path, file_pattern = "sig_genes.RDS"))
#
#   plotted_sig_genes = readRDS(check_for_file(path, file_pattern = "plotSigGenes_data.RDS"))
#
#   # load Z matrix
#   z_mat_path = check_for_file(path, file_pattern = "z_mat.csv")
#
#   if ( is.null(z_mat_path)) {
#     cat("\n Z matrix not found. Check path \n")
#     return()
#   }
#
#   z_mat_full = as.matrix(read.csv(z_mat_path, row.names = 1))
#
#
#   slide_res_path = check_for_file(path, file_pattern = "slide_res.RDS")
#
#   all_LFs = NULL
#
#   if ( !is.null(sig_genes_data) ) {
#     all_LFs = stringr::str_to_upper(unique(readRDS(slide_res_path)$marginal_vars))
#
#   } else {
#     cat("\n Couldn't find significant latent factors. Plotting all latent factors \n")
#     all_LFs = colnames(z_mat_full)
#   }
#
#
#
#   if (F) {
#     #
#     #   z_mat = z_mat_full[, all_LFs]
#     #
#     #   z_mat_nonsig = z_mat_full[, setdiff(colnames(z_mat_full), all_LFs)]
#
#     # z_mat = as.matrix(read.csv(z_mat_path, row.names = 1))[, all_LFs]
#     #
#     #
#     # z_mat_nonsig = as.matrix(read.csv(z_mat_path, row.names = 1))[, -c(all_LFs)]
#
#
#     # Have X, Y and Z matrix loaded. Can plot latent factors now
#
#     # # first attach Y to Z matrix
#     # z_combined = cbind.data.frame(df$y, z_mat)
#     #
#     # names(z_combined)[1] = "y"
#     #
#     # z_nonsig = cbind.data.frame(df$y, z_mat_nonsig)
#     # names(z_nonsig)[1] = "y"
#     #
#     #
#     # cliff_slide = apply(z_combined[, -1], 2, function(x) cliff.delta(f = x, d = z_combined$y)$estimate)
#     # cliff_df = data.frame(cliff_slide)
#     # colnames(cliff_df)[1] = "cliff"
#     # cliff_df$LF = paste0("SZ", 1:length(rownames(cliff_df)))
#     # cliff_df$condition = "SLIDE"
#     #
#     # cliff_nonsig = apply(z_nonsig[, -1], 2, function(x) cliff.delta(f = x, d = z_nonsig$y)$estimate)
#     # cliff_nonsigdf = data.frame(cliff_nonsig)
#     # colnames(cliff_nonsigdf)[1] = "cliff"
#     # cliff_nonsigdf$LF = rownames(cliff_nonsigdf)
#     # cliff_nonsigdf$condition = "Non-sig"
#     #
#     # cliff_combined = rbind(cliff_df, cliff_nonsigdf)
#
#     # z_pivot = z_combined %>% tidyr::pivot_longer(.,
#     #                                              cols = colnames(z_mat),
#     #                                              names_to = "LF_num",
#     #                                              values_to = "LF_val")
#     #
#     # # rename latent factors
#     # for (z in 1:length(all_LFs)) {
#     #
#     #   z_pivot$LF_num = ifelse(z_pivot$LF_num == all_LFs[z], z, z_pivot$LF_num)
#     # }
#     #
#     # z_pivot$LF_val = as.numeric(z_pivot$LF_val)
#     # z_pivot$y = as.numeric(z_pivot$y)
#     #
#     # if (!is.null(group_names) && length(group_names) >= 2) {
#     #   z_pivot$y = ifelse(z_pivot$y == min(z_pivot$y), group_names[1], group_names[2])
#     #
#     #   z_pivot$y = factor(z_pivot$y)
#     #
#     # }
#   }
#
#   dir_name = paste0(path, "/LF_cliffs_delta")
#
#   if ( !dir.exists(dir_name) ) {
#     dir.create(dir_name)
#   }
#   original_wd = getwd()
#
#   if (!is.null(LF_to_plot)) {
#     all_LFs = all_LFs[LF_to_plot]
#   }
#
#   z_mat = z_mat_full[, all_LFs]
#
#   # df$y = df$y[sample(1:length(df$y))]
#   z_combined = cbind.data.frame(df$y, z_mat)
#   names(z_combined)[1] = "y"
#
#
#   # get gene data
#   sig_genes_vec = unlist(lapply(sig_genes_data, function(x) x['gene']))
#
#   sig_genes_plotted_vec = plotted_sig_genes$name
#
#   gene_mat = df[, which(!colnames(df) %in% sig_genes_vec)]
#
#   # sample_size = ifelse(length(sig_genes_vec) > ncol(gene_mat), ncol(gene_mat), length(sig_genes_vec))
#   sample_size = length(sig_genes_plotted_vec)
#
#   sampled_non_sig_genes = sample(1:ncol(gene_mat), sample_size)
#   # sample DF
#   gene_mat = df[, sampled_non_sig_genes]
#
#   gene_mat = cbind.data.frame(df$y, gene_mat)
#   names(gene_mat)[1] = "y"
#
#
#   # run cliffs delta
#
#   doing_classification = ifelse(length(unique(df$y)) == 2, T, F)
#
#   if (doing_classification) {
#
#     z_combined$y = as.numeric(z_combined$y)
#
#     # cliff_slide = abs(apply(z_combined[, -1], 2, function(x) cliff.delta(d = x, f = z_combined$y)$estimate))
#     cliff_slide = abs(apply(z_combined[, -1], 2, function(x) cliff.delta(as.numeric(x[z_combined$y == 0]), as.numeric(x[z_combined$y == 1]))$estimate))
#
#
#     cliff_df = data.frame(cliff_slide)
#     colnames(cliff_df)[1] = "cliff"
#     cliff_df$source = paste0("SZ", 1:length(rownames(cliff_df)))
#     cliff_df$condition = "SLIDE"
#
#     # cliff_genes = abs(apply(gene_mat[, -1], 2, function(x) cliff.delta(d = x, f = gene_mat$y)$estimate))
#     cliff_genes = abs(apply(gene_mat[, -1], 2, function(x) cliff.delta(as.numeric(x[gene_mat$y == 0]), as.numeric(x[gene_mat$y == 1]))$estimate))
#
#     cliff_df2 = data.frame(cliff_genes)
#     colnames(cliff_df2)[1] = "cliff"
#     cliff_df2$source = ""#colnames(gene_mat[, -1])
#     cliff_df2$condition = "LF genes"
#
#     plot_lab = "Cliffs Delta"
#
#   } else {
#     cliff_slide = apply(z_combined[, -1], 2, function(x) abs(cor(x, z_combined$y)))
#     cliff_df = data.frame(cliff_slide)
#     colnames(cliff_df)[1] = "cliff"
#     cliff_df$source = paste0("SZ", 1:length(rownames(cliff_df)))
#     cliff_df$condition = "SLIDE"
#
#     cliff_genes = apply(gene_mat[, -1], 2, function(x) abs(cor(x, gene_mat$y)))
#     cliff_df2 = data.frame(cliff_genes)
#     colnames(cliff_df2)[1] = "cliff"
#     cliff_df2$source = ""#colnames(gene_mat[, -1])
#     cliff_df2$condition = "LF genes"
#
#     plot_lab = "Correlation with Y"
#
#   }
#
#
#
#   cliff_combined = rbind(cliff_df, cliff_df2)
#
#
#   saveRDS(list(cliffs_delta = cliff_combined,
#                sampled_non_sig_genes = sampled_non_sig_genes),
#           paste0(dir_name, "/cliff_results.RDS"))
#
#   # write.csv(sampled_non_sig_genes, paste0(dir_name, "/sampled_non_sig_genes.csv"))
#
#
#   # z_pivot = z_combined %>% tidyr::pivot_longer(.,
#   #                                              cols = colnames(z_mat),
#   #                                              names_to = "LF_num",
#   #                                              values_to = "LF_val")
#
#
#
#   pl = ggpubr::ggboxplot(cliff_combined, x = "condition", y = "cliff",
#                          xlab = "Group", repel = T, fill = "condition",
#                          palette = "aaas",
#                          ylab = plot_lab, add = "point", label = "source") +
#     ggpubr::stat_compare_means(method = "wilcox.test", label.x.npc = "right")
#
#
#   plot_title = paste0(dir_name, "/CD_boxplot_LF_vs_nonsig.png")
#
#   ggplot2::ggsave(filename = plot_title,
#                   plot = pl, height = 7, width = 7,
#                   device = "png")
#
#   #   # for(LF in unique(z_pivot$LF_num)) {
#   #     LF_temp = z_pivot %>% filter(LF_num == LF)
#   #
#   #
#   #     # plots for regression case
#   #     is_classification = length(unique(z_pivot$y)) == 2
#   #
#   #     if (is_classification) {
#   #
#   #
#   #       # if (!is.null(group_names) && length(group_names) >= 2) {
#   #       #   LF_temp$y = ifelse(LF_temp$y == min(LF_temp$y), group_names[1], group_names[2])
#   #       # }
#   #
#   #       pl = ggpubr::ggboxplot(LF_temp,
#   #                              x = "y",
#   #                              y = "LF_val", palette = "npg",
#   #                              xlab = "Group",
#   #                              ylab = paste0("Latent Factor #", LF) ) +
#   #         ggpubr::stat_compare_means(method = "wilcox.test")
#   #
#   #       plot_title = paste0(dir_name, "/LF", LF, "_boxplot.png")
#   #       ggplot2::ggsave(filename = plot_title,
#   #                       plot = pl, height = 7, width = 7,
#   #                       device = "png")
#   #
#   #       # save all LF plots in one
#   #
#   #       pl = ggpubr::ggboxplot(z_pivot,
#   #                              x = "y",
#   #                              y = "LF_val",
#   #                              facet.by = "LF_num",
#   #                              palette = "npg",
#   #                              xlab = "Group",
#   #                              ylab = "Latent Factor Value" ) +
#   #         ggpubr::stat_compare_means(method = "wilcox.test")
#   #
#   #       plot_title = paste0(dir_name, "/", "all_boxplot.png")
#   #       ggplot2::ggsave(filename = plot_title,
#   #                       plot = pl, height = 10, width = 10,
#   #                       device = "png")
#   #
#   #
#   #     } else {
#   #
#   #       # sort by y value
#   #       # LF_temp = LF_temp[order(LF_temp$y), ]
#   #       LF_temp$y = round(as.numeric(LF_temp$y), 3)
#   #
#   #       pl = ggpubr::ggscatter(LF_temp, x = "y",
#   #                              y = "LF_val",
#   #                              palette = "npg",
#   #                              xlab = "Scaled LoSAI",
#   #                              add = "reg.line",
#   #                              add.params = list(color = "blue", fill = "lightgray"),
#   #                              conf.int = T,
#   #                              ylab = paste0("Latent Factor #", LF) ) +
#   #         stat_cor(method = "pearson")
#   #
#   #       plot_title = paste0(dir_name, "/LF", LF, "_regplot.png")
#   #       ggplot2::ggsave(filename = plot_title,
#   #                       plot = pl, height = 7, width = 7,
#   #                       device = "png")
#   #
#   #       pl = ggpubr:: ggscatter(z_pivot, x = "y",
#   #                               y = "LF_val",
#   #                               palette = "npg",
#   #                               xlab = "Scaled LoSAI",
#   #                               facet.by = "LF_num",
#   #                               add = "reg.line",
#   #                               add.params = list(color = "blue", fill = "lightgray"),
#   #                               conf.int = T,
#   #                               ylab = "Significant Latent Factors" ) +
#   #         stat_cor(method = "pearson")
#   #
#   #       plot_title = paste0(dir_name, "/", "all_regplot.png")
#   #       ggplot2::ggsave(filename = plot_title,
#   #                       plot = pl, height = 10, width = 10,
#   #                       device = "png")
#   #
#   #     }
#   #   }
#   #   setwd(original_wd)
#   #
# }
#
#
# adult_ls_vs_healthy = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
#   label = 'Blue = Increased in Healthy \n Red = Increased in LS',
#   group_names = c("Healthy", "LS"), LF_to_plot = c(4, 5) )
#
# all_ls_vs_healthy = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
#   label = 'Blue = Increased in Healthy \n Red = Increased in LS',
#   group_names = c("Healthy", "LS"), LF_to_plot = c(1, 4))
#
#
# ls_losai = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
#   label = 'Blue = Lower LoSAI \n Red = Higher LoSAI',
#   group_names = NULL, LF_to_plot = c(1, 4))
#
# ls_onset = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
#   label = "Blue = Increased in Peds, Red = Increased in Adults",
#   group_names = c("Peds", "Adult"), LF_to_plot = c(4, 5))
#
# for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_losai, ls_onset)) {
#   plot_LF_boxplots(m$path, group_names = m$group_names)
# }













