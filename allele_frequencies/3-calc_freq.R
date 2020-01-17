library(igraph)
library(tidyverse)

all_profiles <- read_tsv("../input_data/integrated_data.tsv")   
    
obs_alleles <- all_profiles %>%
    pivot_longer(m_1:af_2, names_to = "id", values_to = "allele") %>%
    distinct(marker, allele)

obs_alleles_integers <- filter(obs_alleles, near(allele %% 1, 0))
obs_alleles_microvars <- filter(obs_alleles, !near(allele %% 1, 0))

possible_alleles <- obs_alleles_integers %>%
    group_by(marker) %>%
    complete(allele = full_seq(allele, 1)) %>%
    ungroup() %>%
    bind_rows(obs_alleles_microvars) %>%
    arrange(marker, allele)

profilesdf <- all_profiles %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    pivot_longer(m_1:af_2, 
                 names_pattern = "(.+)_(1|2)",
                 names_to = c("person", "h")) %>%
    mutate(h = paste0("h", h)) %>%
    pivot_wider(names_from = h, values_from = value) %>%
    unite(id, c("case_no", "person"), sep = "_")

dups <- sprintf("/scratch/vitor/str/duplicates_%s.tsv", 1:100) %>%
    map_df(read_tsv) %>%
    mutate(id1 = pmin(ind1, id),
           id2 = pmax(ind1, id)) %>%
    distinct(id1, id2)

g <- graph_from_edgelist(as.matrix(dups), directed = FALSE)

dup_groups_df <- tibble(id = names(components(g)$membership), 
                        grp = components(g)$membership)

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
    pivot_longer(-marker, names_to = c("h"), values_to = "allele") %>%
    count(marker, allele) %>%
    group_by(marker) %>%
    mutate(f = n/sum(n)) %>%
    ungroup()

freq_alleles <- select(freq_df, marker, allele)

all_alleles_freq <- left_join(possible_alleles, freq_df, by = c("marker", "allele")) %>%
    mutate(n = replace_na(n, 1)) %>%
    group_by(marker) %>%
    mutate(f = n/sum(n)) %>%
    ungroup()

write_tsv(all_alleles_freq, "./allele_frequencies.tsv")
