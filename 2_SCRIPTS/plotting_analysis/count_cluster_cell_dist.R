library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)

library(ggpubr)

## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')


clust_names = data.frame(table(data@active.ident))

pl = ggpubr::ggbarplot(clust_names, x = "Var1", y = "Freq", fill = "Var1", palette = c(ggpubr::get_palette("aaas", 10), ggpubr::get_palette("npg", 5))) +
  rotate_x_text() + guides(color = guide_legend(nrow = 3))
ggsave(plot = pl, filename = "/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/cell_clust_counts.png",
       width = 8, height = 8)
# %>% as.data.frame() %>% summarise(n = n())
