library(tidyverse)

source('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/3_PLOTTING_SCRIPTS/CalcCliffDelta.R')

library(effsize)
library(tidyverse)
library(ggpubr)

plot_LF_boxplots = function(path, df_filename_RDS_pattern = "df.RDS",
                            group_names = NULL) {


  check_for_file = function(path, file_pattern) {
    f = list.files(path, recursive = F, full.names = T, pattern = file_pattern)
    if( length(f) > 0 ) {
      return(f[1])
    } else {
      cat("\n Folder must have exactly one file with desired name: ", file_pattern, "\n")
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
  z_mat_path = check_for_file(path, file_pattern = "z_mat.csv")

  if ( is.null(z_mat_path)) {
    cat("\n Z matrix not found. Check path \n")
    return()
  }

  z_mat_full = as.matrix(read.csv(z_mat_path, row.names = 1))

  sig_genes_data = check_for_file(path, file_pattern = "slide_res.RDS")


  all_sig_genes = unlist(lapply(readRDS(paste0(path, "/sig_genes.RDS")), function(x) x['gene']))

  gene_idx = which(colnames(df) %in% all_sig_genes)

  all_LFs = NULL

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(unique(readRDS(sig_genes_data)$marginal_vars))

  } else {
    cat("\n Couldn't find significant latent factors. Plotting all latent factors \n")
    all_LFs = colnames(z_mat_full)
  }

  z_mat = z_mat_full[, all_LFs]


  y = as.matrix(df$y)

  z_mat_nonsig = z_mat_full[, setdiff(colnames(z_mat_full), all_LFs)]

  nonsig_inds = 1:length(colnames(z_mat_nonsig))

  # calculate cliffs delta for slide LFs
  lf_inds = as.numeric(stringr::str_split_i(all_LFs, pattern = "Z", i = 2))

  slide_cd = CalcCliffDelta(z_mat, y, comb = list(c(0, 1)), sig_idx = lf_inds)



  nonsig_inds = as.numeric(stringr::str_split_i(colnames(z_mat_nonsig), pattern = "Z", i = 2))

  alt_lf_cd = CalcCliffDelta(z_mat_nonsig, y, comb = list(c(0, 1)), sig_idx = nonsig_inds)


  permute_pvals = CliffDeltaPermute(z_mat_full, y, 50,
                                    slide_cd$all_cds, alt_lf_cd$all_cds,
                                    comb = list(c(0, 1)),
                                    sig_idx = lf_inds)
  # z_mat = z_mat_full[, all_LFs]
  #
  #
  # y = as.matrix(df$y)
  #
  # z_mat_nonsig = z_mat_full[, setdiff(colnames(z_mat_full), all_LFs)]
  #
  # nonsig_inds = 1:length(colnames(z_mat_nonsig))
  #
  # # calculate cliffs delta for slide LFs
  # lf_inds = as.numeric(stringr::str_split_i(all_LFs, pattern = "Z", i = 2))
  #
  # slide_cd = CalcCliffDelta(z_mat, y, comb = list(c(0, 1)), sig_idx = lf_inds)
  #
  #
  #
  # nonsig_inds = as.numeric(stringr::str_split_i(colnames(z_mat_nonsig), pattern = "Z", i = 2))
  #
  # alt_lf_cd = CalcCliffDelta(z_mat_nonsig, y, comb = list(c(0, 1)), sig_idx = nonsig_inds)
  #
  #
  # permute_pvals = CliffDeltaPermute(z_mat_full, y, 50,
  #                                   slide_cd$all_cds, alt_lf_cd$all_cds,
  #                                   comb = list(c(0, 1)),
  #                                   sig_idx = lf_inds)



  all_results = list(slide_cd = slide_cd,
                     alt_lf_cd = alt_lf_cd,
                     permute_pvals = permute_pvals)


  # aggregate the dfs

  slide_df = slide_cd$all_cds
  alt_lf_df = alt_lf_cd$all_cds
  null_df = permute_pvals[[3]]$full_null_cds


  slide_df$method = "SLIDE"
  alt_lf_df$method = "Nonsig LF"
  null_df$method = "Null"

  dir_name = paste0(path, "/LF_cliffs_delta")

  full_df = rbind(rbind(slide_df, alt_lf_df), null_df)


  pl = ggpubr::ggboxplot(full_df, x = "method", y = "deltas", fill = "method",
                         palette = "aaas",
                         xlab = "Latent Factor", repel = T,
                         ylab = "Cliffs Delta", add = "point") +
    ggpubr::stat_compare_means(method = "wilcox.test", ref.group = "SLIDE")

  plot_title = paste0(dir_name, "/CD_boxplot_HX.png")

  ggplot2::ggsave(filename = plot_title,
                  plot = pl, height = 7, width = 7,
                  device = "png")



  # z_mat_nonsig = z_mat_full[, setdiff(colnames(z_mat_full), all_LFs)]
  #
  # # Have X, Y and Z matrix loaded. Can plot latent factors now
  #
  # # first attach Y to Z matrix
  # z_combined = cbind.data.frame(df$y, z_mat)
  #
  # names(z_combined)[1] = "y"
#
#   z_nonsig = cbind.data.frame(df$y, z_mat_nonsig)
#   names(z_nonsig)[1] = "y"
#
#
#   cliff_slide = apply(z_combined[, -1], 2, function(x) cliff.delta(f = x, d = z_combined$y)$estimate)
#   cliff_df = data.frame(cliff_slide)
#   colnames(cliff_df)[1] = "cliff"
#   cliff_df$LF = paste0("SZ", 1:length(rownames(cliff_df)))
#   cliff_df$condition = "SLIDE"
#
#   cliff_nonsig = apply(z_nonsig[, -1], 2, function(x) cliff.delta(f = x, d = z_nonsig$y)$estimate)
#   cliff_nonsigdf = data.frame(cliff_nonsig)
#   colnames(cliff_nonsigdf)[1] = "cliff"
#   cliff_nonsigdf$LF = rownames(cliff_nonsigdf)
#   cliff_nonsigdf$condition = "Non-sig"

  # cliff_combined = rbind(cliff_df, cliff_n onsigdf)

  # dir_name = paste0(path, "/LF_cliffs_delta")
  #
  # if ( !dir.exists(dir_name) ) {
  #   dir.create(dir_name)
  # }
  # original_wd = getwd()
  #
  # pl = ggpubr::ggboxplot(cliff_combined, x = "condition", y = "cliff",
  #                        xlab = "Latent Factor", repel = T,
  #                        ylab = "Cliffs Delta", add = "point", label = "LF") +
  #   ggpubr::stat_compare_means(method = "wilcox.test")
  #
  # plot_title = paste0(dir_name, "/CD_boxplot.png")
  #
  # ggplot2::ggsave(filename = plot_title,
  #                 plot = pl, height = 7, width = 7,
  #                 device = "png")

  saveRDS(all_results, paste0(path, "/cliffs_delta_results.RDS"))

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

for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset)) {
  plot_LF_boxplots(m$path, group_names = m$group_names)
}


