library(tidyverse)

freqs <- read_tsv("../input_data/allele_frequency.tsv")

simul_expr <- expression(freqs %>% 
    group_by(marker) %>% 
    sample_n(size = 2, replace = TRUE, weight = f) %>% 
    ungroup() %>%
    select(-f))

set.seed(1)
simul_af <- replicate(1e6, eval(simul_expr), simplify = FALSE) %>%
    bind_rows(.id = "case_no")
  
simul_ch <- replicate(1e6, eval(simul_expr), simplify = FALSE) %>%
    bind_rows(.id = "case_no")

simul_df <- list(af = simul_af, ch = simul_ch) %>%
    bind_rows(.id = "subject") %>%
    arrange(case_no, marker, subject, allele) %>%
    group_by(case_no, subject, marker) %>%
    mutate(h = 1:2) %>%
    ungroup() %>%
    unite(id, c("subject", "h"), sep = "_") %>%
    spread(id, allele)

write_tsv(simul_df, "./simulated_duos.tsv")