## this code makes correlation networks ##
# load packages
library(tidyverse)
library(qgraph) ## for making the network
# library(patchwork)



# igv_palette = c("#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338",
#   "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#837B8D", "#D58F5C",
#   "#7A65A5", "#E4AF69", "#3B1B53")
#
#
# # # #   "chr18" = "#CDDEB7",
  # "chr19" = "#612A79", "chr20" = "#AE1F63", "chr21" = "#E7C76F",
  # "chr22" = "#5A655E", "chrX" = "#CC9900", "chrY" = "#99CC00",
  # "chrUn" = "#A9A9A9", "chr23" = "#CC9900", "chr24" = "#99CC00",
  # "chr25" = "#33CC00", "chr26" = "#00CC33", "chr27" = "#00CC99",
  # "chr28" = "#0099CC", "chr29" = "#0A47FF", "chr30" = "#4775FF",
  # "chr31" = "#FFC20A", "chr32" = "#FFD147", "chr33" = "#990033",
  # "chr34" = "#991A00", "chr35" = "#996600", "chr36" = "#809900",
  # "chr37" = "#339900", "chr38" = "#00991A", "chr39" = "#009966",
  # "chr40" = "#008099", "chr41" = "#003399", "chr42" = "#1A0099",
  # "chr43" = "#660099", "chr44" = "#990080", "chr45" = "#D60047",
  # "chr46" = "#FF1463", "chr47" = "#00D68F", "chr48" = "#14FFB1"
# # # )

make_color_map = function(cluster_ids_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/condensed_clusters.txt') {
  data = read_delim(cluster_ids_path, delim = " ")

  # append palette
  # igv_palette = c("#5050FF", "#CE3D32", "#749B58", "#466983", "#F0E685", "#BA6338",
  #                 "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#E7C76F", "#D58F5C",
  #                 "#D60047", "#E4AF69", "#CBC3E3")#"#660099") #"#3B1B53")
  igv_palette = c("#5050FF", "#CE3D32", "#749B58", "#466983", "#F0E685", "#BA6338",
                  "#5DB1DD", "#802268", "#003399", "#D595A7", "#E7C76F", "#D58F5C",
                  "#D60047", "#E4AF69", "#CBC3E3")#"#660099") #"#3B1B53")
  data$color = igv_palette[1]

  all_cells = unique(data$Celltype)

  for (u in 1:length(all_cells)) {

    data[data$Celltype == all_cells[u], "color"] = igv_palette[u]
  }

  return(data)
}

get_color_from_cluster = function(genes) {


  # cluster_data = read_delim('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/prev_hardcoded_clusters.txt', delim = " ")

  cluster_data = read_delim('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/hardcoded_clusters_new_palette.txt', delim = " ")

  # get cluster nums from gene names
  clust_nums = stringr::str_extract(genes, pattern = "\\.[0-9]+\\.")

  # remove period
  clust_nums = as.numeric(stringr::str_remove_all(clust_nums, pattern = "\\."))


  # data_order = match(cluster_data$Seurat_Clusters, clust_nums)
  # clust_colors = cluster_data[cluster_data$Seurat_Clusters %in% clust_nums, "color"]
  clust_colors = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "color"])#cluster_data[cluster_data$Seurat_Clusters %in% clust_nums, "color"]

  # clust_colors = unlist(clust_colors)

  clust_names = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "Celltype"])

  return(list(clust_names = clust_names, clust_colors = clust_colors))
}


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

  dir_name = paste0(path, "/cell_type_corr_network")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }

  original_wd = getwd()
  pl_list = list()
  for ( LF in unique(sig_genes_data$fac) ) {
    sg_temp = ( sig_genes_data %>% filter(fac == LF) )

    clust_res = get_color_from_cluster(sg_temp$names)
    clust_names = unlist(clust_res$clust_names)
    clust_colors = unlist(clust_res$clust_colors)


    clust_names_list = list()
    color_names_list = list()
    color_names_vec = c()
    for (a in unique(clust_names)) {

      clust_names_list[[a]] = unlist(which(clust_names == a))

      color_names_vec = c(color_names_vec, unique(unlist(clust_colors[which(clust_names == a)])))

    }

    x_cor = x[, sg_temp$names]

    # remove C.number from gene names
    colnames(x_cor) = stringr::str_split_i(colnames(x_cor), pattern = "C\\.[0-9]+\\.", i = 2)

    x_temp = cor(as.matrix(x_cor))

    node_shape = ifelse(is.na(sg_temp$color), "circle", ifelse(sg_temp$color == "Red", "triangle", "square"))

    out_dir = paste0(dir_name, "/LF", LF, "_corr")

    # change wd cause weird
    setwd(dir_name)
    plot_title = ifelse(is.null(plot_label), "", plot_label)

    # plot_title = paste0("LF #", LF, "\n", plot_title, "\nTriangle = Up\nSquare = Down")
    plot_title = paste0("LF #", LF)

    # pl = qgraph(x_temp, filename = paste0("celltype_LF_shape_cluster", LF), groups = clust_names_list,
    #             color = color_names_vec,
    #             shape = node_shape,
    #             layout = "groups", threshold=0.5, repulsion=0.1,
    #             labels = colnames(x_temp), GLratio = 3.5,
    #             # title = plot_title,
    #             label.scale.equal=TRUE,label.prop=0.95,shape="ellipse",
    #             posCol="#51087E", negCol="#59A14F",filetype='pdf',height=5,width=8)
    # make layout using shapes
    pl_shape = qgraph(x_temp, filename = paste0("celltype_LF_shape_cluster", LF), groups = node_shape,
                nodeNames = clust_names,
                legend.mode = "groups",
                color = clust_colors,
                shape = node_shape,
                layout = "spring",
                # minimum=0.1,
                # repulsion=1.5,
                labels = colnames(x_temp),
                # GLratio = 3.5,
                # title = plot_title,
                label.scale.equal=TRUE,label.prop=0.8,shape="ellipse",
                posCol="#51087E", negCol="#59A14F",filetype='pdf',height=10,width=10)

    # pl_shape$layout
    pl = qgraph(x_temp, filename = paste0("celltype_LF_shape_cluster_fade_", LF), groups = clust_names_list,
                # legend.mode = "groups",
                color = color_names_vec,
                shape = node_shape,
                edge.width = 1,
                node.width = 0.75,
                node.height = 0.75, fade = T,
                # layout = "groups",
                # threshold=0.3,
                minimum=0.25,
                # repulsion=1,
                labels = colnames(x_temp),
                # GLratio = 3.5,
                # title = plot_title,
                layout = pl_shape$layout,

                label.scale.equal=TRUE,label.prop=0.95,shape="ellipse",
                posCol="#51087E", negCol="#59A14F",filetype='pdf',height=5,width=8)

    pl_list[[LF]] = pl
  }
  saveRDS(pl_list, paste0(dir_name, "/plot_list.RDS"))

  setwd(original_wd)
}



adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  # label = 'Blue = Increased in Healthy \nRed = Increased in LS')
  label = "")

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  # label = 'Blue = Increased in Healthy \nRed = Increased in LS')
  label = "")

ls_losai = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
  # label = 'Blue = Lower LoSAI \nRed = Higher LoSAI')
  label = "")

ls_onset = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
  label = "")


# for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_losai, ls_onset)) {
for (m in list(all_ls_vs_healthy)) {

  create_corr_network(m$path, m$label)
}

