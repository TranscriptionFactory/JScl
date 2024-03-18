library(tidyverse)
library(ggpubr)
metadata <- readRDS("/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/data/metadata.RDS")


######################
# sex_plot =
# sex_health_plot =


health_onset = ggplot(metadata, aes(x = health, fill = onset)) + geom_bar(position = "dodge") +
  scale_fill_manual(values = ggpubr::get_palette("aaas", 2)) + theme_pubr()
ggsave(filename = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/onset_health.png',
       plot = health_onset, width = 6, height = 5)

sex_health = ggplot(metadata, aes(x = health, fill = sex)) + geom_bar(position = "dodge") +
  scale_fill_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) + theme_pubr()
ggsave(filename = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/sex_health.png',
plot = sex_health, width = 6, height = 5)


subtype_sex = ggplot(metadata %>% filter(health == "LS"), aes(x = subtype, fill = sex)) + geom_bar(position = "dodge") +
    scale_fill_manual(values = ggpubr::get_palette("aaas", 4)[3:4]) + ggtitle("LS Subtype") + theme_pubr(x.text.angle = 90)
ggsave(filename = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/subtype_sex.png',
plot = subtype_plot, width = 8, height = 7)

subtype_onset = ggplot(metadata %>% filter(health == "LS"), aes(x = subtype, fill = onset)) + geom_bar(position = "dodge") +
  scale_fill_manual(values = ggpubr::get_palette("aaas", 4)[1:2]) + ggtitle("LS Subtype") + theme_pubr(x.text.angle = 90)
ggsave(filename = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/subtype_onset.png',
       plot = subtype_onset, width = 8, height = 7)



subtype_losai = ggplot(metadata %>% filter(health == "LS") %>% mutate(LoSAI.mLoSSI = as.numeric(LoSAI.mLoSSI)) %>%
                         mutate(PGA.A = as.numeric(PGA.A)), aes(x = LoSAI.mLoSSI, y = PGA.A, color = factor(subtype))) + geom_point(size = 5) +
  scale_color_manual(values = ggpubr::get_palette("aaas", 20)) + ggtitle("LS Subtype") + theme_pubr(x.text.angle = 90)
ggsave(filename = '/ix/djishnu/Aaron/3_collabs/1_Torok_Jscl_ER/ER_analysis/5_manuscript_figures/metadata_plots/subtype_LoSAI.png',
       plot = subtype_losai, width = 8, height = 7)
