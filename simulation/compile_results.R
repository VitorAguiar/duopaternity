library(tidyverse)


#r001 

dat_r001_all <- read_tsv("./results/cpi_r001_all.tsv")
dat_r001_codis <- read_tsv("./results/cpi_r001_codis.tsv")
dat_r001_ident <- read_tsv("./results/cpi_r001_ident.tsv")
dat_r001_pp16 <- read_tsv("./results/cpi_r001_pp16.tsv")

inclusion_r001_codis <- dat_r001_codis %>%
    summarise(n = mean(n_exclusions < 3 & cpi >= 10000)*100)

inclusion_r001_ident <- dat_r001_ident %>%
    summarise(n = mean(n_exclusions < 3 & cpi >= 10000)*100)

inclusion_r001_pp16 <- dat_r001_pp16 %>%
    summarise(n = mean(n_exclusions < 3 & cpi >= 10000)*100)


#STRBASE

dat_strbase_all <- read_tsv("./results/cpi_strbase_all.tsv")
dat_strbase_ident <- read_tsv("./results/cpi_strbase_ident.tsv")
dat_strbase_pp16 <- read_tsv("./results/cpi_strbase_pp16.tsv")

inclusion_strbase_ident <- dat_strbase_ident %>%
    summarise(n = mean(n_exclusions < 3 & cpi >= 10000)*100)

inclusion_strbase_pp16 <- dat_strbase_pp16 %>%
    summarise(n = mean(n_exclusions < 3 & cpi >= 10000)*100)

