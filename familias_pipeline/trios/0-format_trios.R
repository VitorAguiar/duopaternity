library(tidyverse)

freq <- read_tsv("../../allele_frequencies/allele_frequencies.tsv")

trios <- read_tsv("../../input_data/integrated_data.tsv") %>%
    filter(marker %in% freq$marker) %>%
    group_by(case_no, trio) %>%
    filter(n() > 18) %>%
    ungroup()

trios_familias_format <- trios %>%
    pivot_longer(m_1:af_2, names_to = c("person", "h"), names_pattern = "(.+)_(.)", values_to = "allele") %>%
    mutate(person = recode(person, "m" = "mother", "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    pivot_wider(names_from = m, values_from = allele)

trios_batch <- trios_familias_format %>%
    distinct(case_no, trio) %>%
    mutate(batch = ntile(1:nrow(.), 8))

trios_out <- trios_familias_format %>%
    left_join(trios_batch, by = c("case_no", "trio")) %>%
    select(case_no, trio, batch, everything())

write_tsv(trios_out, "./trios_familias_format.tsv")
