library(tidyverse)

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_exc <- read_tsv("../trios/trios_cpi_eq.tsv") %>%
    filter(n_loci >= 18, n_exclusions >= 4) %>%
    select(case_no, trio) %>%
    inner_join(trios, by = c("case_no", "trio"))

duos <- trios_exc %>%
    select(-m_1, -m_2)  

duos_familias_format <- duos %>%
    pivot_longer(ch_1:af_2, 
                 names_to = c("person", "h"), 
                 names_pattern = "(.+)_(.)", 
                 values_to = "allele") %>%
    mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    pivot_wider(names_from = m, values_from = allele)

write_tsv(duos_familias_format, "duos_familias_format.tsv")
