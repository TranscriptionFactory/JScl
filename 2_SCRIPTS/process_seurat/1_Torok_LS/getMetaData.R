
library(tidyverse)

metadata = read.csv('/ix/djishnu/Aaron/Jscl_ER/scRNAseq_LS_Update_Metadata_Jacob.csv')

# look at severity within Peds LS
y_data = read.csv('/ix/djishnu/Aaron/Jscl_ER/INPUT_DATA/check_yPedsLS.csv')

# get corresponding column for PGA

pedsLS_rows = which(metadata$sample_id %in% y_data$library_id)

yj = y_data %>% left_join(metadata, by = c("library_id" = "sample_id"))

new_y = yj %>% select(library_id, "PGA.A")
names(new_y)[2] = "health"

write.csv(new_y, '/ix/djishnu/Aaron/Jscl_ER/INPUT_DATA/check_yPedsLS_pga.csv',
          row.names = F)


# look at severity within adult LS
y_data = read.csv('/ix/djishnu/Aaron/Jscl_ER/INPUT_DATA/check_yAdultLS.csv')

# get corresponding column for PGA

adultLS_rows = which(metadata$sample_id %in% y_data$library_id)

yj = y_data %>% left_join(metadata, by = c("library_id" = "sample_id"))

new_y = yj %>% select(library_id, "PGA.A")
names(new_y)[2] = "health"

write.csv(new_y, '/ix/djishnu/Aaron/Jscl_ER/INPUT_DATA/check_yAdultLS_pga.csv',
          row.names = F)
