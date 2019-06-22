library(tidyverse)

total_trios <- read_tsv("../familias_pipeline/trios/total_trios.tsv") %>%
    select(1, 3) %>%
    gather(i, n, -1) %>%
    spread(marker_set, n) %>%
    select(-i)

#r001 

dat_r001_all <- read_tsv("./results/cpi_r001_all.tsv")
dat_r001_codis <- read_tsv("./results/cpi_r001_codis.tsv")
dat_r001_codisplus <- read_tsv("./results/cpi_r001_codisplus.tsv")
dat_r001_ident <- read_tsv("./results/cpi_r001_ident.tsv")
dat_r001_pp16 <- read_tsv("./results/cpi_r001_pp16.tsv")

inclusion_r001_codis <- dat_r001_codis %>%
    slice(1:(total_trios$codis*1000)) %>%
    mutate(slot = rep(1:total_trios$codis, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_codis_r001 <- inclusion_r001_codis %>% 
    count(slot) %>%
    filter(n >= obs_r001$codis) %>%
    count()

inclusion_r001_codisplus <- dat_r001_codisplus %>%
    slice(1:(total_trios$codisplus*1000)) %>%
    mutate(slot = rep(1:total_trios$codisplus, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_codis_r001 <- inclusion_r001_codis %>% 
    count(slot) %>%
    filter(n >= obs_r001$codis) %>%
    count()

inclusion_r001_ident <- dat_r001_ident %>%
    slice(1:(total_trios$identifiler*1000)) %>%
    mutate(slot = rep(1:total_trios$identifiler, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_ident_r001 <- inclusion_r001_ident %>% 
    count(slot) %>%
    filter(n >= obs_r001$identifiler) %>%
    count()

inclusion_r001_pp16 <- dat_r001_pp16 %>%
    slice(1:(total_trios$pp16*1000)) %>%
    mutate(slot = rep(1:total_trios$pp16, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_pp16_r001 <- inclusion_r001_pp16 %>% 
    count(slot) %>%
    filter(n >= obs_r001$pp16) %>%
    count()

#STRBASE

dat_strbase_all <- read_tsv("./results/cpi_strbase_all.tsv")
dat_strbase_ident <- read_tsv("./results/cpi_strbase_ident.tsv")
dat_strbase_pp16 <- read_tsv("./results/cpi_strbase_pp16.tsv")

inclusion_strbase_ident <- dat_strbase_ident %>%
    slice(1:(total_trios$identifiler*1000)) %>%
    mutate(slot = rep(1:total_trios$identifiler, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_ident_strbase <- inclusion_strbase_ident %>% 
    count(slot) %>%
    filter(n >= obs_strbase$identifiler) %>%
    count()

inclusion_strbase_pp16 <- dat_strbase_pp16 %>%
    slice(1:(total_trios$pp16*1000)) %>%
    mutate(slot = rep(1:total_trios$pp16, each = 1000)) %>%
    group_by(slot) %>%
    filter(n_exclusions < 3, cpi >= 10000) %>%
    ungroup()

simul_pp16_strbase <- inclusion_strbase_pp16 %>% 
    count(slot) %>%
    filter(n >= obs_strbase$pp16) %>%
    count()

