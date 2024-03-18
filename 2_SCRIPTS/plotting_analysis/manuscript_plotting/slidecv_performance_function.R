require(doParallel)
require(dplyr)
library(SLIDE)
library(EssReg)
library(yaml)
library(ggpubr)

args <- commandArgs(trailingOnly = T)
run_path  <- args[1]
#replicate  <- args[2]

# get yaml
yaml_params = readRDS(paste0(run_path, "/yaml_parameters.RDS"))

yaml_params$x_path = paste0(run_path, "/x.csv")
yaml_params$y_path = paste0(run_path, "/y.csv")
yaml_params$z_path = paste0(run_path, "/z_mat.csv")
yaml_params$parallel = TRUE
yaml_params$fdr = 0.1
yaml_params$f_size = 20
yaml_params$spec = 0.1
yaml_params$niter = 1000

yaml_path = paste0(run_path, "/yaml_params.yaml")

# make new dir for slide results
slide_output_dir = paste0(run_path, "/slide_cv/")
if (!dir.exists(slide_output_dir)) {
  dir.create(slide_output_dir, recursive = T)
}

yaml_params$out_path = slide_output_dir

yaml::write_yaml(yaml_params, yaml_path)

er_input <- yaml::yaml.load_file(yaml_path)

for(j in 1:er_input$nreps){

  withCallingHandlers(SLIDE::paperBench(yaml_path, replicate=j))

}


pathLists <- list.files(er_input$out_path,recursive = T,pattern = "results")
perfList <- lapply(paste0(er_input$out_path,pathLists), readRDS)
perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


if (er_input$eval_type == "corr") {
  lambda_boxplot <- ggplot2::ggplot(data = perRes,
                                    ggplot2::aes(x = method,
                                                 y = corr,
                                                 fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(fill = "Method") +
    ggplot2::scale_alpha(guide = 'none')
} else {
  lambda_boxplot <- ggplot2::ggplot(data = perRes,
                                    ggplot2::aes(x = method,
                                                 y = auc,
                                                 fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(fill = "Method") +
    ggplot2::scale_alpha(guide = 'none')
}
library(ggplot2)
ggsave(paste0(er_input$out_path,"/delta",er_input$delta,"lambda",er_input$lambda,"_boxplot.pdf"),lambda_boxplot)
saveRDS(perRes,file=paste0(er_input$out_path,"/delta",er_input$delta,"lambda",er_input$lambda,"_boxplot_data.rds"))
