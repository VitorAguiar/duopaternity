library(tidyverse)

codis_loci <- readLines("../input_data/codis_loci.txt")
codisplus_loci <- readLines("../input_data/codisplus_loci.txt")
ident_loci <- readLines("../input_data/identifiler_loci.txt")
pp16_loci <- readLines("../input_data/pp16_loci.txt")

simul_r001 <- paste0("./results/pi_r001_chunk", 1:25, ".tsv") %>%
    map_df(read_tsv)

cpi_r001_all <- simul_r001 %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_r001_all, "./results/cpi_r001_all.tsv")

cpi_r001_codis <- simul_r001 %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_r001_codis, "./results/cpi_r001_codis.tsv")

cpi_r001_codisplus <- simul_r001 %>%
    filter(marker %in% codisplus_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_r001_codisplus, "./results/cpi_r001_codisplus.tsv")

cpi_r001_ident <- simul_r001 %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_r001_ident, "./results/cpi_r001_ident.tsv")

cpi_r001_pp16 <- simul_r001 %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi), 
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_r001_pp16, "./results/cpi_r001_pp16.tsv")

