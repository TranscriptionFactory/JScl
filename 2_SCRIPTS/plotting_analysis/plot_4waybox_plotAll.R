# peds signal in all
library(tidyverse)


## this code makes correlation networks ##
# load packages
library(tidyverse)
library(ggpubr) ## for making the network
# library(patchwork)



# igv_palette = c("#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338",
#   "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#837B8D", "#D58F5C",
#   "#7A65A5", "#E4AF69", "#3B1B53")
#
#
# # # #   "chr18" = "#CDDEB7",
#   "chr19" = "#612A79", "chr20" = "#AE1F63", "chr21" = "#E7C76F",
#   "chr22" = "#5A655E", "chrX" = "#CC9900", "chrY" = "#99CC00",
#   "chrUn" = "#A9A9A9", "chr23" = "#CC9900", "chr24" = "#99CC00",
#   "chr25" = "#33CC00", "chr26" = "#00CC33", "chr27" = "#00CC99",
#   "chr28" = "#0099CC", "chr29" = "#0A47FF", "chr30" = "#4775FF",
#   "chr31" = "#FFC20A", "chr32" = "#FFD147", "chr33" = "#990033",
#   "chr34" = "#991A00", "chr35" = "#996600", "chr36" = "#809900",
#   "chr37" = "#339900", "chr38" = "#00991A", "chr39" = "#009966",
#   "chr40" = "#008099", "chr41" = "#003399", "chr42" = "#1A0099",
#   "chr43" = "#660099", "chr44" = "#990080", "chr45" = "#D60047",
#   "chr46" = "#FF1463", "chr47" = "#00D68F", "chr48" = "#14FFB1"
# # # )




get_peds_vs_adult = function(data,
                             metadata_path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/metadata.RDS'){
  metadata = readRDS(metadata_path)

  metadata = metadata[metadata$sample_id %in% data, c("sample_id", "onset", "health")]
  rownames(metadata) = metadata$sample_id
  return(metadata)
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

get_clust_info = function(genes) {

  cluster_data = read_delim('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/prev_hardcoded_clusters.txt', delim = " ")

  # get cluster nums from gene names
  clust_nums = stringr::str_extract(genes, pattern = "\\.[0-9]+\\.")

  gene_name_no_clust = stringr::str_split_i(genes, pattern = "C\\.[0-9]+\\.", i = 2)

  # remove period
  clust_nums = as.numeric(stringr::str_remove_all(clust_nums, pattern = "\\."))

  clust_colors = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "color"])#cluster_data[cluster_data$Seurat_Clusters %in% clust_nums, "color"]

  clust_names = sapply(clust_nums, function(x) cluster_data[cluster_data$Seurat_Clusters == as.numeric(x), "Celltype"])




  df = cbind.data.frame(unlist(clust_names), gene_name_no_clust, genes)
  colnames(df) = c("clust", "genes", "x_col")
  return(list(clust_names = clust_names, clust_colors = clust_colors, genes = genes, gene_name = gene_name_no_clust,
              df = df))
}


create_corr_network = function(path, encodings, comp = NULL, comp_var = NULL, plot_label = NULL) {

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


  # sig_genes_data = check_for_file(path, file_pattern = "plotSigGenes_data.RDS")
  sig_genes_data = check_for_file(path, file_pattern = "sig_genes.RDS")

  sig_genes_data = readRDS(sig_genes_data)

  all_LFs = c()

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(sig_genes_data$names)

  } else {
    cat("\n Couldn't find genes in significant latent factors. Check path. \n")
    return()
  }

  # color_code = stringr::str_to_lower(sig_genes_data$color)

  x = df[, -1]

  if ( dim(x)[2] <= 1) {
    cat("\n Latent factor genes not found in X data. Check path \n")
    return()
  }

  dir_name = paste0(path, "/keratinocyte")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }

  original_wd = getwd()

  all_sg = unlist(lapply(sig_genes_data, function(x) x['gene']))
  # all_sg = sig_genes_data$names


  # indx = stringr::str_detect(all_sg, pattern = "SKP1|COX7C|EEF1B2|ACTG1", negate = F)
  #
  # if (length(which(indx)) != 0) {
  #   all_sg = all_sg[stringr::str_detect(all_sg, pattern = "SKP1|COX7C|EEF1B2|ACTG1", negate = T)]
  #
  # }


  sg_df = get_clust_info(all_sg)$df

  # select these genes from x
  x_keratinocyte = df[, c("y", sg_df$x_col)]

  x_keratinocyte$y = ifelse(x_keratinocyte$y == 0, encodings[1], encodings[2])

  x_keratinocyte[, comp_var] = get_peds_vs_adult(rownames(x))[, comp_var]

  x_keratinocyte$y = factor(x_keratinocyte$y, levels = encodings)

  x_long = x_keratinocyte %>% pivot_longer(cols = sg_df$x_col, names_to = "gene")


  x_long$clust = sapply(x_long$gene, function(x) sg_df[sg_df$x_col == x, "clust"][1])
  x_long$gene = sapply(x_long$gene, function(x) sg_df[sg_df$x_col == x, "genes"][1])




  # x_subset = x_long[stringr::str_detect(x_long$clust, "Keratinocytes"), ]
  x_subset = x_long#[stringr::str_detect(x_long$clust, "Keratinocytes"), ]

  for(c in unique(x_subset$clust)) {

    temp = x_subset %>% filter(clust == c)


    # peds variance
    temp_peds = temp %>% filter(comp_var == comp) %>% group_by(gene) %>% summarize(gene_var = var(value))

    temp_peds = temp_peds[order(temp_peds$gene_var, decreasing = T), ]

    # take top 4
    temp = temp %>% filter(gene %in% unique(temp_peds$gene)[1:4])

    # p = ggpubr::ggboxplot(temp, x = "y", y = "value", fill = comp_var, facet.by = "gene", palette = "aaas",
    #                       ylab = "Expression", xlab = "Group", title = c, ncol = 4, outlier.shape = NA,
    #                       scales = "free_y") + rotate_x_text() + stat_compare_means(comparisons = list(c("Healthy", "LS"), c("Peds", "Adult")))

    p = ggpubr::ggboxplot(temp, x = "y", y = "value", fill = comp_var, facet.by = "gene", palette = "aaas",
                          ylab = "Expression", xlab = "Group", title = c, ncol = 4, outlier.shape = NA,
                          scales = "free_y") + rotate_x_text() + stat_compare_means(comparisons = list(c("Healthy", "LS")))
    plot_length = ifelse(length(unique(temp_peds$gene)) < 4, 2.5 * length(unique(temp_peds$gene)), 10)
    ggplot2::ggsave(filename = paste0(dir_name, "/", c, "_boxplot_all_clust_no_outlier.pdf"), plot = p, width = plot_length, height = 4.25)

    # plot_height = ceiling(length(unique(temp$gene))/4)
    # plot_height = ifelse(plot_height == 1, (5/3), plot_height)
    # ggplot2::ggsave(filename = paste0(dir_name, "/", c, "_boxplot_filt.pdf"), plot = p, width = 10, height = 3  * plot_height)

  }


  # p = ggpubr::ggboxplot(x_subset, x = "y", y = "value", fill = "onset", facet.by = "gene", palette = "aaas",
  #                       ylab = "Expression", xlab = "Group",
  #                       scales = "free_y") + rotate_x_text()
  #
  # # ggplot2::ggsave(filename = paste0(dir_name, "/4way_boxplot.png"), plot = p, width = 10, height = 10) + stat_compare_means()
  # ggplot2::ggsave(filename = paste0(dir_name, "/4way_boxplot.pdf"), plot = p, width = 10, height = 10)

}



adult_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/Adult_LS_vs_Healthy',
  # label = 'Blue = Increased in Healthy \nRed = Increased in LS')
  encodings = c("Healthy", "LS"),
  comp = "LS",
  comp_var = "health")

all_ls_vs_healthy = list(
  path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy',
  # label = 'Blue = Increased in Healthy \nRed = Increased in LS')
  encodings = c("Healthy", "LS"),
  comp = "Peds",
  comp_var = "onset")


# ls_losai = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_LoSAI',
#   # label = 'Blue = Lower LoSAI \nRed = Higher LoSAI')
#   label = "")

# ls_onset = list(
#   path = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/LS_onset',
#   encodings = c("Peds", "Adult"))
#


# for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy, ls_onset)) {
for (m in list(adult_ls_vs_healthy, all_ls_vs_healthy)) {

  create_corr_network(m$path, m$encodings, m$comp, m$comp_var)
}

