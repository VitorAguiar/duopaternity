library(tidyverse)
library(cowplot)
library(scales)


trios_cpi <- read_tsv("./familias_pipeline/trios/trios_cpi_sw.tsv") %>%
    select(-n_loci, -n_exclusions)
    
trios_exclusions <- read_tsv("./familias_pipeline/trios/trios_cpi_eq.tsv") %>%
    select(case_no, trio, n_exclusions)

trios_df <- left_join(trios_cpi, trios_exclusions, by = c("case_no", "trio"))
duos_df <- read_tsv("./familias_pipeline/duos/duos_cpi.tsv")


cpi_df <- left_join(duos_df, trios_df, 
                    by = c("case_no", "trio"),
                    suffix = c(".Duos", ".Trios")) %>%
    pivot_longer(-(1:2),
                 names_pattern = "(.+)\\.(.+)",
                 names_to = c("stat", "case_type")) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(case_type = factor(case_type, levels = c("Trios", "Duos")))
    

plot1 <- ggplot(cpi_df, aes(case_type, cpi)) +
    geom_hline(yintercept = 1e5) +
    geom_hline(yintercept = 1e4, linetype = 2) +
    geom_blank(aes(y = 1e10)) +
    geom_boxplot() +
    scale_y_log10("CPI",
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(family = "Times", size = 12),
          axis.text = element_text(family = "Times", size = 11))

plot2 <- ggplot(cpi_df, aes(case_type, n_exclusions)) +
    geom_boxplot() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(family = "Times", size = 12),
          axis.text = element_text(family = "Times", size = 11)) +
    labs(y = "# of excluding loci")
    

plot_out <- plot_grid(plot1, plot2, nrow = 1, 
                      labels = c("A)", "B)"), 
                      label_size = 14, label_fontfamily = "Times", label_x = -0.025)

save_plot("./Fig1.tiff", plot_out, dpi = 600, base_width = 6, base_height = 2.5)
