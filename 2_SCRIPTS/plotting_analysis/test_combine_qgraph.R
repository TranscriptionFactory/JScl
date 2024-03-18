library(qgraph)
library(tidyverse)
plot_list <- readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/2_FINAL_RESULTS/All_LS_vs_Healthy/cell_type_corr_network/plot_list.RDS")

# plot.new()
par(mfrow = c(2,2))

nf = layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE), respect = TRUE)

for (p in 1:length(plot_list)) {
  tplot = ifelse(p == 1, TRUE, FALSE)
  qgraph::qgraph(p, plot = T)
}
