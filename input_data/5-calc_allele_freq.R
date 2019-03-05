library(tidyverse)

profiles <- read_tsv("./profiles.tsv", guess_max = 2e5)

res <- readRDS("./dnatools_results.rds")

resdf <- as_tibble(res$hits) %>%
    mutate_at(vars(id1, id2), as.character) %>% 
    mutate(total_match = match + Fmatch) %>%
    filter(total_match >= 21) %>%
    arrange(total_match) %>%
    mutate(pair = 1:n()) %>%
    select(pair, id1, id2)

to_del <- resdf %>%
    gather(ix, id, id1:id2) %>%
    left_join(profiles) %>%
    gather(marker, allele, -(pair:id)) %>%
    separate(marker, c("marker", "h"), sep = "\\.") %>%
    arrange(pair, marker, ix, h) %>%
    group_by(pair, id) %>%
    mutate(nas = sum(allele == 0)) %>%
    ungroup() %>%
    select(pair:marker, nas) %>%
    distinct(pair, ix, nas, .keep_all = TRUE) %>%
    group_by(pair) %>%
    slice(which.max(nas)) %>%
    ungroup() %>%
    distinct(id) %>%
    pull(id)

profiles_noDups <- profiles %>%
    filter(! id %in% to_del) %>%
    gather(marker, allele, -id) %>%
    filter(allele != 0) %>%
    separate(marker, c("marker", "h"), sep = "\\.")

allele_counts <- profiles_noDups %>%
    select(marker, allele) %>% 
    add_count(marker) %>%
    add_count(marker, allele, name = "nn") %>%
    distinct()

mafs <- allele_counts %>%
    distinct(marker, n) %>%
    mutate(maf = 5/n) %>%
    select(marker, maf)

parents_freqs <- allele_counts %>%
    mutate(f = nn/n) %>%
    select(marker, allele, f)

alleles_not_in_parents <- read_tsv("./integrated_data.tsv") %>%
    gather(h, allele, -(1:3)) %>%
    select(marker, allele) %>%
    anti_join(parents_freqs) %>%
    arrange(marker, allele) %>%
    mutate(f = 0)

freqs <- bind_rows(parents_freqs, alleles_not_in_parents) %>%
    left_join(mafs) %>%
    mutate(f = ifelse(f < maf, maf, f)) %>%
    group_by(marker) %>%
    mutate(f = f/sum(f)) %>%
    ungroup() %>%
    select(marker, allele, f) %>%
    arrange(marker, allele)

write_tsv(freqs, "./allele_frequency.tsv")
