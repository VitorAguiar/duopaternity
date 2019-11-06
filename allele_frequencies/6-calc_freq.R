library(tidyverse)

microvars <- read_tsv("../input_data/str_info.tsv") %>%
    mutate(mv = map(repeats, ~0:(.x-1))) %>%
    unnest(mv) %>%
    select(-repeats)

all_profiles <- read_tsv("../input_data/integrated_data.tsv")   

possible_alleles <- all_profiles %>%
    pivot_longer(m_1:af_2, names_to = "id", values_to = "allele") %>%
    distinct(marker, allele) %>%
    group_by(marker) %>%
    summarise(minallele = min(allele) - 1,
	      maxallele = max(allele) + 1) %>%
    ungroup() %>%
    mutate(minallele = ifelse(minallele <= 1, 1, minallele),
	   minallele = floor(minallele),
	   maxallele = floor(maxallele)) %>%
    mutate(allelerange = map2(minallele, maxallele, ~.x:.y)) %>%
    select(marker, allelerange) %>%
    unnest(allelerange) %>%
    left_join(microvars, by = c("marker" = "locus")) %>%
    mutate(allele = ifelse(mv == 0, allelerange, paste(allelerange, mv, sep = ".")),
	   allele = as.numeric(allele)) %>%
    select(marker, allele)

profilesdf <- all_profiles %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    gather(id, allele, 3:6) %>%
    separate(id, c("ind", "h"), sep = "_") %>%
    unite(id, c("case_no", "ind"), sep = "_") %>%
    mutate(h = paste0("h", h)) %>%
    spread(h, allele)

dups <- read_tsv("./matchtbl_all.tsv") %>%
    mutate(id01 = pmin(id1, id2),
	   id02 = pmax(id1, id2)) %>%
    distinct(id01, id02) %>%
    rename(id1 = id01, id2 = id02)

dup_groups <- vector("list", nrow(dups))

for (i in 1:nrow(dups)) {
    
    ids <- slice(dups, i) %>% unlist()

    dup_groups[[i]] <- dups %>% 
	filter(id1 %in% ids | id2 %in% ids) %>%
	unlist() %>%
	sort() %>%
	unique()
}

dup_groups_df <- unique(dup_groups) %>%
    map(~tibble(id = .)) %>%
    bind_rows(.id = "grp")

# calculate allele frequencies
dups_to_calc <- dup_groups_df %>%
    left_join(profilesdf) %>%
    distinct(grp, marker, h1, h2) %>% 
    add_count(grp, marker) %>% 
    filter(n == 1) %>%
    select(-grp, -n)

non_dups_to_calc <- profilesdf %>%
    filter(! id %in% dup_groups_df$id) %>%
    select(-id)

freq_df <- bind_rows(non_dups_to_calc, dups_to_calc) %>%
    gather(h, allele, h1:h2) %>%
    count(marker, allele) %>%
    group_by(marker) %>%
    mutate(f = n/sum(n)) %>%
    ungroup()

freq_alleles <- select(freq_df, marker, allele)

all_alleles_freq <- left_join(possible_alleles, freq_df, by = c("marker", "allele")) %>%
    mutate(n = replace_na(n, 2)) %>%
    group_by(marker) %>%
    mutate(f = n/sum(n)) %>%
    ungroup()

write_tsv(all_alleles_freq, "./allele_frequencies.tsv")
