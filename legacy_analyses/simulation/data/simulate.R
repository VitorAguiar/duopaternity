library(tidyverse)

inclusion <- read_tsv("../../familias_pipeline/duos/duos_inclusion_stepwise.tsv")

duos_sw <- read_tsv("../../familias_pipeline/duos/duos_pi_stewise.tsv") %>%
    filter(case_no %in% inclusion$case_no)

freqs <- read_tsv("../../allele_frequencies/allele_frequencies.tsv") %>%
    select(marker, allele, f = adj_f) %>%
    filter(marker %in% duos_sw$marker)

n_cases <- 1e6

set.seed(1)

simul_duos <- freqs %>% 
	group_by(marker) %>% 
	sample_n(size = n_cases * 4, replace = TRUE, weight = f) %>%
	mutate(case_no = rep(1:n_cases, each = 4)) %>%
	ungroup() %>%
	mutate(ind = rep(c("af", "af", "ch", "ch"), n()/4L)) %>%
	arrange(case_no, ind, marker, allele) %>%
	mutate(h = rep(1:2, n()/2L)) %>%
	select(case_no, marker, ind, h, allele) %>% 
	unite(id, c("ind", "h"), sep = "_") %>%
	spread(id, allele) 

write_tsv(simul_duos, "./simulated_duos.tsv")
