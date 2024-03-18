library(tidyverse)
library(effsize)

# x, y are input data frames
# returns list of genes sorted by cliffs delta (high to low)
filter_by_cliffs_deta = function(x, y) {

  x = as.matrix(x)
  y = as.matrix(y)

  cliff_res = apply(x, 2, function(x_temp) abs(cliff.delta(f = x_temp, d = y)$estimate))

  cliff_res = data.frame("cliff_delta" = cliff_res)
  cliff_res$name = rownames(cliff_res)

  # sort
  cliff_res = cliff_res[order(cliff_res$cliff_delta, decreasing = T), ]
  return(cliff_res)
}


make_color_map = function(cluster_ids_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/condensed_clusters.txt') {
  data = read_delim(cluster_ids_path, delim = " ")

  # append palette
  igv_palette = c("#5050FF", "#CE3D32", "#749B58", "#466983", "#F0E685", "#BA6338",
                  "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#E7C76F", "#D58F5C",
                  "#D60047", "#E4AF69", "#CBC3E3")#"#660099") #"#3B1B53")

  data$color = igv_palette[1]

  all_cells = unique(data$Celltype)

  for (u in 1:length(all_cells)) {

    data[data$Celltype == all_cells[u], "color"] = igv_palette[u]
  }

  return(data)
}


get_color_from_cluster = function(genes) {

  # cluster_data = make_color_map()
  cluster_data = read_delim('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/hardcoded_clusters_new_palette.txt', delim = " ")

  # get cluster nums from gene names
  clust_nums = stringr::str_extract(genes, pattern = "\\.[0-9]+\\.")

  gene_name = stringr::str_split_i(genes, pattern = "\\.[0-9]+\\.", i = 2)
  # remove period
  clust_nums = as.numeric(stringr::str_remove_all(clust_nums, pattern = "\\."))

  # data_order = match(cluster_data$Seurat_Clusters, clust_nums)
  # clust_colors = cluster_data[cluster_data$Seurat_Clusters %in% clust_nums, "color"]
  clust_colors = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "color"])#cluster_data[cluster_data$Seurat_Clusters %in% clust_nums, "color"]
  clust_colors = unlist(clust_colors)

  clust_names = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "Celltype"])
  return(list(clust_names = clust_names, clust_colors = clust_colors, gene_name = gene_name))
}



get_LF_data = function(er_run_path, lfs) {

  er_data = readRDS(list.files(er_run_path, pattern = "plotSigGenes_data.RDS", full.names = T)[1])
  # select LFs

  er_data_subset = er_data[er_data$fac %in% unique(er_data$fac)[lfs], ]

  return(er_data_subset)
  # get effect sizes
}


get_x_and_y = function(path) {
  check_for_file = function(path, file_pattern) {

    f = list.files(path, recursive = T, full.names = T, pattern = file_pattern)
    if( length(f) > 0 ) {
      return(f[1])
    } else {
      cat("\n Folder must have exactly one file with desired name: ", pattern, "\n")
      return(NULL)
    }
  }


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


subset_big_df = function(df, metadata, lfs) {


  df = df[, colnames(df) %in% lfs]

  # get effect size ranking

  return(df)
}



make_bubbleplot = function(data, sig_gene_data, metadata, out_path = NULL) {


  if (!is.null(out_path)) {
    out_path_folder = paste0(out_path, "/bubble_plots/")
    if (!dir.exists(out_path_folder)) {
      dir.create(out_path_folder)

    }
    out_path = out_path_folder

  }

  results = list()

  for (LF in unique(sig_gene_data$fac)) {

    genes = sig_gene_data[sig_gene_data$fac == LF, ]

    # pivot data


    found_genes = base::intersect(colnames(data), genes$names)

    metahealth = as.matrix(metadata[,"health"])
    metahealth = ifelse(metahealth == "LS", 1, 0)
    found_genes = filter_by_cliffs_deta(data[, found_genes], metahealth)$name[1:10]
    df_wide = cbind(data[, found_genes], metadata)
    df_long = df_wide %>% tidyr::pivot_longer(cols = found_genes, names_to = "gene",
                                              values_to = "exp")

    df_long$group = paste0(df_long$onset, "_", df_long$health)

    # add average expression column
    df_long = df_long %>% group_by(group, gene) %>% summarise(avg_exp = mean(exp))

    clust_name_colors = get_color_from_cluster(df_long$gene)

    df_long$gene_name = clust_name_colors$gene_name
    df_long$cluster = as.factor(unlist(clust_name_colors$clust_names))
    df_long$clust_colors = clust_name_colors$clust_colors


    p = ggplot(df_long, aes(x = gene_name, y = group, size = avg_exp, color = cluster)) +
      geom_count(show.legend = c("size" = F, "color" = T)) +
      xlab("") + ylab("") +
      # scale_color_manual(values = unique(df_long$clust_colors)) +
      scale_color_manual(values = setNames(df_long$clust_colors,  df_long$cluster)) +

      scale_size(range=c(2,10)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90), legend.position = "top",  legend.text = element_text(size = 12)) +
      guides(color = guide_legend(nrow = 3)) + #theme(text = element_text(size = 12)) +
      ggtitle(paste0("LF ", LF))

    results[[length(results) + 1]] = p

    if (!is.null(out_path)) {
      ggsave(filename = paste0(out_path, "LF_new_big_bubble_", LF, ".pdf"), p, height = 5, width = 7)
    }
  }
  return(results)
}




make_bubbleplot_all_lsvsh = function(data, sig_gene_data, metadata, out_path = NULL, cliff_filter = 0, cliff_y = NULL) {


  if (!is.null(out_path)) {
    out_path_folder = paste0(out_path, "/bubble_plots/")
    if (!dir.exists(out_path_folder)) {
      dir.create(out_path_folder)

    }
    out_path = out_path_folder

  }

  results = list()

  for (LF in unique(sig_gene_data$fac)) {

    genes = sig_gene_data[sig_gene_data$fac == LF, ]

    # pivot data

    found_genes = base::intersect(colnames(data), genes$names)
    # found_genes = genes$names
    # if (cliff_filter > 0) {
    #
    # metahealth = as.matrix(metadata[,"health"])
    # metahealth = ifelse(metahealth == "LS", 1, 0)
    # found_genes = filter_by_cliffs_deta(data[, found_genes], metahealth)$name[1:cliff_filter]
    # }
    # uncomment for All LS vs H
    metahealth = as.matrix(metadata[,"health"])
    metahealth = ifelse(metahealth == "LS", 1, 0)
    # found_genes = filter_by_cliffs_deta(data[, found_genes], metahealth)$name[1:12]
    found_genes = filter_by_cliffs_deta(data[, found_genes], metahealth)$name


    df_wide = cbind(data[, found_genes], metadata)
    df_long = df_wide %>% tidyr::pivot_longer(cols = found_genes, names_to = "gene",
                                              values_to = "exp")

    df_long$group = paste0(df_long$onset, "_", df_long$health)

    # add average expression column
    df_long = df_long %>% group_by(group, gene) %>% summarise(avg_exp = mean(exp))

    clust_name_colors = get_color_from_cluster(df_long$gene)

    df_long$gene_name = clust_name_colors$gene_name
    df_long$cluster = as.factor(unlist(clust_name_colors$clust_names))
    df_long$clust_colors = clust_name_colors$clust_colors


    p = ggplot(df_long, aes(x = gene_name, y = group, size = avg_exp, color = cluster)) +
      geom_count(show.legend = c("size" = F, "color" = T)) +
      xlab("") + ylab("") +
      scale_color_manual(values = setNames(df_long$clust_colors,  df_long$cluster)) +
      scale_size(range=c(2,10)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90), legend.position = "right",  legend.text = element_text(size = 12)) +
      # guides(color = guide_legend(nrow = 4)) + #theme(text = element_text(size = 12)) +
      ggtitle(paste0("LF ", LF))

    results[[length(results) + 1]] = p

    if (!is.null(out_path)) {
      ggsave(filename = paste0(out_path, "LF_new_big_bubble_all_right_legend", LF, ".pdf"), p, height = 5, width = 7)
    }
  }
  return(results)
}



make_bubbleplot_lsonset = function(data, sig_gene_data, metadata, out_path = NULL) {


  if (!is.null(out_path)) {
    out_path_folder = paste0(out_path, "/bubble_plots/")
    if (!dir.exists(out_path_folder)) {
      dir.create(out_path_folder)

    }
    out_path = out_path_folder

  }

  results = list()

  for (LF in unique(sig_gene_data$fac)) {

    genes = sig_gene_data[sig_gene_data$fac == LF, ]

    # pivot data


    found_genes = base::intersect(colnames(data), genes$names)

    # metahealth = as.matrix(metadata[,"onset"])
    # metahealth = ifelse(metahealth[,1] == "Adult", 1, 0)
    # found_genes = filter_by_cliffs_deta(data[, found_genes], metahealth)$name
    # if (length(found_genes) > 10) {
    #   found_genes = found_genes[1:10]
    # }
    # found_genes = ifelse(length(found_genes) > 10, found_genes[1:10], found_genes[1:length(found_genes)])

    df_wide = cbind(data[, found_genes], metadata)
    df_long = df_wide %>% tidyr::pivot_longer(cols = found_genes, names_to = "gene",
                                              values_to = "exp")

    df_long$group = paste0(df_long$onset, "_", df_long$health)

    # add average expression column
    df_long = df_long %>% group_by(group, gene) %>% summarise(avg_exp = mean(exp))

    clust_name_colors = get_color_from_cluster(df_long$gene)

    df_long$gene_name = clust_name_colors$gene_name
    df_long$cluster = as.factor(unlist(clust_name_colors$clust_names))
    df_long$clust_colors = clust_name_colors$clust_colors


    p = ggplot(df_long, aes(x = gene_name, y = group, size = avg_exp, color = cluster)) +
      geom_count(show.legend = c("size" = F, "color" = T)) +
      xlab("") + ylab("") +
      # scale_color_manual(values = unique(df_long$clust_colors)) +
      scale_color_manual(values = setNames(df_long$clust_colors,  df_long$cluster)) +

      scale_size(range=c(2,10)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90), legend.position = "right",  legend.text = element_text(size = 12)) +
      # guides(color = guide_legend(nrow = 3)) + #theme(text = element_text(size = 12)) +
      ggtitle(paste0("LF ", LF))

    results[[length(results) + 1]] = p

    if (!is.null(out_path)) {
      ggsave(filename = paste0(out_path, "LF_new_big_bubble_all_right_legend", LF, ".pdf"), p, height = 5, width = 7)
    }
  }
  return(results)
}

