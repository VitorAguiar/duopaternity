library(tidyverse)

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_cpi <- read_tsv("./trios_cpi.tsv") %>%
    filter(n_exclusions < 4, cpi > 100000) %>%
    select(case_no, trio)

trios_m <- inner_join(trios_cpi, trios) %>%
    select(case_no, trio, marker, ch_1, ch_2, af_1, af_2) %>%
    group_by(marker) %>%
    summarise(m = mean(af_1 != ch_1 & af_1 != ch_2 & af_2 != ch_1 & af_2 != ch_2))

aabb <- read_tsv("../../input_data/aabb_mutation_rates.tsv")



left_join(trios_m, aabb) %>% print(n = Inf)

inner_join(trios_cpi, trios) %>%
    select(case_no, trio, marker, ch_1, ch_2, af_1, af_2) %>%
    filter(af_1 != ch_1 & af_1 != ch_2 & af_2 != ch_1 & af_2 != ch_2)
