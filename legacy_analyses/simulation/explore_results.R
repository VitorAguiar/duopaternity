library(tidyverse)

duos <- data.table::fread("./data/simulated_duos.tsv") %>%
    as_tibble()

freqs <- read_tsv("../allele_frequencies/allele_frequencies.tsv")

pi_r001 <- data.table::fread("./duos_pi_r001.tsv") %>% 
    as_tibble()

pi_sw <- data.table::fread("./duos_pi_stepwise.tsv") %>% 
    as_tibble()

left_join(pi_r001, pi_sw, by = c("case_no", "marker"), suffix = c(".r001", ".sw")) %>%
    select(-starts_with("exclusion")) %>%
    filter(case_no == inc_r001$case_no[3]) %>% print(n = Inf)

duos %>% filter(case_no == inc_r001$case_no[1]) %>% print(n = Inf)

freqs %>% filter(marker == "D22S1045") %>% print(n = Inf)

inc_r001 <- read_tsv("./duos_inclusion_r001.tsv")
inc_sw <- read_tsv("./duos_inclusion_stepwise.tsv")


full_join(inc_r001, inc_sw, by = "case_no", suffix = c(".r001", ".sw")) %>%
    select(case_no, cpi.r001, cpi.sw) %>%
    as.data.frame()

full_join(inc_r001, inc_sw, by = "case_no", suffix = c(".r001", ".sw")) %>%
    mutate(diff_cpi = cpi.r001 - cpi.sw,
           perc_increase = diff_cpi/cpi.sw * 100) 



