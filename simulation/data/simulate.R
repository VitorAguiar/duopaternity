library(tidyverse)

freqs <- read_tsv("../../input_data/allele_frequency.tsv")

n_cases <- 1e6

set.seed(1)

for (i in 1:25) {

    case_no_end <- i * n_cases
    case_no_start <- case_no_end - n_cases + 1

    out <- paste0("./simul_chunk", i, ".tsv") 
    
    freqs %>% 
	group_by(marker) %>% 
	sample_n(size = n_cases * 4, replace = TRUE, weight = f) %>%
	mutate(case_no = rep(case_no_start:case_no_end, each = 4)) %>%
	ungroup() %>%
	mutate(ind = rep(c("af", "af", "ch", "ch"), n()/4)) %>%
	select(-f) %>%
	arrange(case_no, ind, marker, allele) %>%
	mutate(h = rep(1:2, n()/2L)) %>%
	select(case_no, marker, ind, h, allele) %>% 
	unite(id, c("ind", "h"), sep = "_") %>%
	spread(id, allele) %>%
	write_tsv(out)
}
