library(tidyverse)

profilesdf <- read_tsv("../input_data/integrated_data.tsv") %>%
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

# Alleles in children not observed in parents
freq_alleles <- select(freq_df, marker, allele)

all_profiles <- read_tsv("../input_data/integrated_data.tsv")

all_alleles <- all_profiles %>%
    gather(ind, allele, m_1:af_2)

not_in_parents <- anti_join(all_alleles, freq_alleles) %>% distinct(marker, allele)

# Update allele freqs
# alleles with f > 5/2N or alleles not observed in parents will get f = 5/2N
updated_freq_df <- bind_rows(freq_df, not_in_parents) %>%
    group_by(marker) %>%
    mutate(f2 = ifelse(n < 5 | is.na(n), 5/sum(n, na.rm = TRUE), f),
	   adj_f = f2/sum(f2)) %>%
    ungroup() %>%
    arrange(marker, allele)

write_tsv(updated_freq_df, "./allele_frequencies.tsv")
