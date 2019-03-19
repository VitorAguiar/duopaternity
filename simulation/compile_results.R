library(tidyverse)

total_trios <- read_tsv("../familias_pipeline/trios/total_trios.tsv") %>%
    select(1, 3) %>%
    gather(i, n, -1) %>%
    spread(marker_set, n) %>%
    select(-i)

dat_r001_all <- read_tsv("./results/cpi_r001_all.tsv")
dat_r001_codis <- read_tsv("./results/cpi_r001_codis.tsv")
dat_r001_ident <- read_tsv("./results/cpi_r001_ident.tsv")
dat_r001_pp16 <- read_tsv("./results/cpi_r001_pp16.tsv")

inclusion_r001_codis <- dat_r001_codis %>%
    slice(1:(total_trios$codis*1000)) %>%
    mutate(slot = rep(1:total_trios$codis, each = 1000)) %>%
    group_by(slot) %>%
    filter(cpi >= 10000) %>%
    ungroup()

inclusion_r001_ident <- dat_r001_ident %>%
    slice(1:(total_trios$identifiler*1000)) %>%
    mutate(slot = rep(1:total_trios$identifiler, each = 1000)) %>%
    group_by(slot) %>%
    filter(cpi >= 10000) %>%
    ungroup()

inclusion_r001_pp16 <- dat_r001_pp16 %>%
    slice(1:(total_trios$pp16*1000)) %>%
    mutate(slot = rep(1:total_trios$pp16, each = 1000)) %>%
    group_by(slot) %>%
    filter(cpi >= 10000) %>%
    ungroup()



