library(tidyverse)

ident_loci <- readLines("../input_data/identifiler_loci.txt")
pp16_loci <- readLines("../input_data/pp16_loci.txt")

simul_strbase <- paste0("./results/pi_strbase_chunk", 1:25, ".tsv") %>%
    map_df(read_tsv)

cpi_strbase_all <- simul_strbase %>%
    group_by(case_no) %>%
    summarise(cpi = prod(pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_strbase_all, "./results/cpi_strbase_all.tsv")

cpi_strbase_ident <- simul_strbase %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_strbase_ident, "./results/cpi_strbase_ident.tsv")

cpi_strbase_pp16 <- simul_strbase %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no) %>%
    summarise(cpi = prod(pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(cpi_strbase_pp16, "./results/cpi_strbase_pp16.tsv")

