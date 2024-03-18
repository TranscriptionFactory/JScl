library(tidyverse)
library(RColorBrewer)
path = list.dirs(recursive = F)[stringr::str_which(list.dirs(recursive = F), pattern = "new_cluster")]


cluster_defs = read.csv('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/new_cluster_defs.txt', sep = "\t")
cluster_defs$Celltype = stringr::str_replace(cluster_defs$Celltype, pattern = "-| ", replacement = "_")


cluster_defs$cluster_cols = colorRampPalette(colors = c("red", "blue", "yellow"))(length(cluster_defs$Cluster))

all_runs = grep(x = list.dirs(path),
                pattern = "lambd(a)?([0-9]+)(\\.)?([0-9]*)_delt(a)?([0-9]+)(\\.)?([0-9]*)$",
                value = T, ignore.case = T)


runs_with_filename = grep(x = list.files(all_runs, full.names = TRUE,recursive = TRUE),
                          pattern = "plotSigGenes_data.RDS$",
                          value = T, ignore.case = T)

valid_runs = stringr::str_split_i(runs_with_filename, i = 1,
                                  pattern = "plotSigGenes_data.RDS$")

# plot with loading
for (i in valid_runs) {

  cdf = readRDS(paste0(i, "/plotSigGenes_data.RDS"))

  clust_nums = as.numeric(stringr::str_split_i(cdf$names, pattern = "\\.", i = 2))

  cluster_nums = unlist(lapply(clust_nums, function(x) cluster_defs[cluster_defs$Cluster == x, "Celltype"]))

  cdf$cluster_cols = unlist(lapply(clust_nums, function(x) cluster_defs[cluster_defs$Cluster == x, "cluster_cols"]))
  cdf$cluster_nums = cluster_nums
  # cdf$cluster_nums = cluster_defs[as.numeric(stringr::str_split_i(cdf$names, pattern = "\\.", i = 2)), ]$Celltype

  cdf$cluster_gene_names = stringr::str_split_i(cdf$names, pattern = "\\.", i = 3)

  # cdf$loading_anno = ifelse(cdf$A_loading == 1, "*", " ")

  cdf$reg_type = ifelse(cdf$color == "Red", "↑ ", "↓ ")

  cdf$names_anno = paste0(cdf$reg_type, cdf$cluster_gene_names)

  plt = cdf %>% ggplot2::ggplot(., aes(x = factor(fac), y = heights, label = names_anno)) +
    ggplot2::geom_text(aes(color = cluster_nums), size = 5) +
    # ggplot2::scale_color_manual(values = cluster_defs$cluster_cols, guide = "none") +
    theme_void() +
    ggplot2::theme(axis.text.x = element_text(), axis.title.x = element_text(),
                   axis.title.y = element_text(angle = 90)) +
    ggplot2::xlab("Latent Factor") + ylab("Genes Associated with Latent Factor Cluster") + ggplot2::ylim(0, 20) +
    ggplot2::ggtitle("Genes Associated with Latent Factor")


  plot_width = ifelse((length(unique(cdf$fac)) * 4 ) > 14, 15, length(unique(cdf$fac)) * 4 )

  ggplot2::ggsave(plot = plt, device = "png",
                  filename = paste0(i, '/plotSigGenes_cluster.png'), width = plot_width, height = 7)

  plt = cdf %>% ggplot2::ggplot(., aes(x = factor(fac), y = heights, label = names)) +
    ggplot2::geom_text(aes(color = factor(color)), size = 5) +
    ggplot2::scale_color_manual(values = c("blue", "red"), guide = "none") + theme_void() +
    ggplot2::theme(axis.text.x = element_text(), axis.title.x = element_text(),
                   axis.title.y = element_text(angle = 90)) +
    ggplot2::xlab("Latent Factor") + ylab("Genes Associated with Latent Factor Cluster") + ggplot2::ylim(0, 20) +
    ggplot2::ggtitle("Genes Associated with Latent Factor")
  #
  #
  # plot_width = ifelse((length(unique(cdf$fac)) * 4 ) > 14, 15, length(unique(cdf$fac)) * 4 )
  #
  ggplot2::ggsave(plot = plt, device = "png",
                  filename = paste0(i, '/plotSigGenes.png'), width = plot_width, height = 7)

}

# for (i in valid_runs) {
#   cdf = readRDS(paste0(i, "/plotSigGenes_data.RDS"))
#   plt = cdf %>% ggplot2::ggplot(., aes(x = factor(fac), y = heights, label = names)) +
#     ggplot2::geom_text(aes(color = factor(color)), size = 5) +
#     ggplot2::scale_color_manual(values = c("blue", "red"), guide = "none") + theme_void() +
#     ggplot2::theme(axis.text.x = element_text(), axis.title.x = element_text(),
#                    axis.title.y = element_text(angle = 90)) +
#     ggplot2::xlab("Latent Factor") + ylab("Genes Associated with Latent Factor Cluster") + ggplot2::ylim(0, 20) +
#     ggplot2::ggtitle("Genes Associated with Latent Factor")
#
#
#   plot_width = ifelse((length(unique(cdf$fac)) * 4 ) > 15, 15, length(unique(cdf$fac)) * 4 )
#
#   ggplot2::ggsave(plot = plt, device = "pdf",
#                   filename = paste0(i, '/big_plotSigGenes.pdf'), width = plot_width, height = 7)
#
# }
