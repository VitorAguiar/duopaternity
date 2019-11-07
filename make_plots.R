library(tidyverse)

ident_cpi_r001 <- read_tsv("./familias_pipeline/duos/duos_ident_cpi_r001.tsv")
pp16_cpi_r001 <- read_tsv("./familias_pipeline/duos/duos_pp16_cpi_r001.tsv")
codis_cpi_r001 <- read_tsv("./familias_pipeline/duos/duos_codis_cpi_r001.tsv")
codisplus_cpi_r001 <- read_tsv("./familias_pipeline/duos/duos_codisplus_cpi_r001.tsv")

ident_cpi_strb <- read_tsv("./familias_pipeline/duos/duos_ident_cpi_strbase.tsv")
pp16_cpi_strb <- read_tsv("./familias_pipeline/duos/duos_pp16_cpi_strbase.tsv")

cpi_df <- list(r001_15A = ident_cpi_r001,
               r001_15B = pp16_cpi_r001,
               r001_codis = codis_cpi_r001,
               r001_codisplus = codisplus_cpi_r001,
               strb_15A = ident_cpi_strb,
               strb_15B = pp16_cpi_strb) %>%
    bind_rows(.id = "set") %>%
    filter(cpi < 10000)

ggplot(cpi_df, aes(cpi, n_exclusions)) +
    geom_point(aes(color = set))

ggplot(cpi_df, aes(set, cpi)) +
    geom_boxplot()

ggplot(cpi_df, aes(cpi, fill = set)) +
    geom_density()

summary(cpi_df$cpi)
