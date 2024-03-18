## this code makes correlation networks ##
# load packages
library(tidyverse)
library(qgraph) ## for making the network
# library(patchwork)

create_corr_network = function(path, plot_label = NULL) {

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
      return(cbind.data.frame(load_matrix_csv(y), load_matrix_csv(x)))
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


  sig_genes_data = check_for_file(path, file_pattern = "plotSigGenes_data.RDS")
  sig_genes_data = readRDS(sig_genes_data)

  all_LFs = c()
  # sig_genes_data = c()

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(sig_genes_data$names)

  } else {
    cat("\n Couldn't find genes in significant latent factors. Check path. \n")
    return()
  }

  color_code = stringr::str_to_lower(sig_genes_data$color)

  x = df[, -1]

  # df[['color_code']] = color_code
  # x = x[, c("color_code", all_LFs)]

  if ( dim(x)[2] <= 1) {
    cat("\n Latent factor genes not found in X data. Check path \n")
    return()
  }

  dir_name = paste0(path, "/corr_network")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }

  original_wd = getwd()

  for ( LF in unique(sig_genes_data$fac) ) {
    sg_temp = ( sig_genes_data %>% filter(fac == LF) )

    x_cor = x[, sg_temp$names]

    x_temp = cor(as.matrix(x_cor))

    out_dir = paste0(dir_name, "/LF", LF, "_corr")

    # change wd cause weird
    setwd(dir_name)
    plot_title = ifelse(is.null(plot_label), "", plot_label)

    plot_title = paste0("LF #", LF, "\n", plot_title)

    pl = qgraph(x_temp, filename = paste0("LF_thresh", LF), color = color_code,
                layout = "spring", threshold=0.0, repulsion=0.4,
                labels = colnames(x_temp),label.font=2,
                title = plot_title, threshold = 0.25,
                label.scale.equal=TRUE,label.prop=0.95,shape="ellipse",
                posCol="#51087E", negCol="#59A14F",filetype='png',height=5,width=10)
  }

  setwd(original_wd)
}



adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \nRed = Increased in LS')

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  label = 'Blue = Increased in Healthy \nRed = Increased in LS')

ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  label = 'Blue = Lower LoSAI \nRed = Higher LoSAI')

ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "Blue = Increased in Peds, Red = Increased in Adults")


for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_losai, ls_onset)) {
  create_corr_network(m$path, m$label)
}

