library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(patchwork)


## LOAD IN DATA


# load("/ix/djishnu/Aaron/Jscl_ER/LS_all_Update_2023.Rdata")
data = LoadH5Seurat('/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/202306_relabeled_clusters/LS.all.H5Seurat')







make_feature_plots = function(path, save_lab = T) {

  get_sig_genes = function(path) {
    return(readRDS(list.files(path, pattern = "^plotSigGenes_data.RDS$", full.names = T)))
  }

  sigdf = get_sig_genes(path)

  for (lf in unique(sigdf$fac)) {
    snames = (sigdf %>% filter(fac == lf))$names

    snames = stringr::str_split_i(snames, "\\.", i = 3)

    pl = FeaturePlot(data, features = snames, label = save_lab, order = T, combine = T,
                     label.size = 4,
                     # cols = c("lightgray", "red", "blue"),
                     keep.scale = "all")

    if (!dir.exists(paste0(path, "/UMAPS"))) {
      dir.create(paste0(path, "/UMAPS"))
    }


    pl_height = ifelse(length(pl) > 16, 30, 20)
    if (save_lab) {

      ggsave(paste0(path, '/UMAPS/lf_', lf, '.pdf'), wrap_plots(pl), width = 30, height = pl_height,
             device = "pdf")

      ggsave(paste0(path, '/UMAPS/lf_', lf, '.png'), wrap_plots(pl), width = 30, height = pl_height,
             device = "png")

    } else {
      ggsave(paste0(path, '/UMAPS/no_lab_lf_', lf, '.pdf'), wrap_plots(pl), width = 30, height = pl_height,
             device = "pdf")

      ggsave(paste0(path, '/UMAPS/no_lab_lf_', lf, '.png'), wrap_plots(pl), width = 30, height = pl_height,
             device = "png")
    }
  }
}

# Adult LS v Healthy
alsh = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/202206_P50_RESULTS/Adult_LS_vs_Healthy'

make_feature_plots(alsh, save_lab = T)

all = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/202206_P50_RESULTS/All_LS_vs_Healthy'
make_feature_plots(all, save_lab = T)

lsl = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/202206_P50_RESULTS/LS_LoSAI'
make_feature_plots(lsl, save_lab = T)

