library(tidyverse)

strtbl <- read_tsv("./str_parents.tsv")

str_loci_counts <- count(strtbl, id) %>%
    mutate(n = n/2L)

dupstbl <- read_tsv("./dupstable.tsv")

left_join(dupstbl, str_loci_counts, by = c("id.x" = "id")) %>%
    rename(n.x = n) %>%
    left_join(str_loci_counts, by = c("id.y" = "id")) %>%
    rename(n.y = n)


