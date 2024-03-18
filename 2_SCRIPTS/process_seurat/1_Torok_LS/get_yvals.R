library(tidyverse)

get_yvals = function(xpath, metapath = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/scRNAseq_LS_Update_Metadata_Jacob.csv'){
  ymat = data.frame("sample_id" = row.names(read.csv(xpath, row.names = 1)))
  metadata = read.csv(metapath)
  
  return(ymat %>% left_join(metadata, by = c("sample_id")))
}
