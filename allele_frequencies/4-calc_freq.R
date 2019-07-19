library(tidyverse)

profilesdf <- read_tsv("../input_data/integrated_data.tsv") %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    gather(id, allele, 3:6) %>%
    separate(id, c("ind", "h"), sep = "_") %>%
    unite(id, c("case_no", "ind"), sep = "_") %>%
    mutate(h = paste0("h", h)) %>%
    spread(h, allele)

# method below potentially fails when A = B and B = C, but 
# A != C because the overlap of loci is low
dups <- read_tsv("./sharetable_all.tsv") %>%
    select(1:2) %>%
    filter(! id.x %in% id.y) %>%
    mutate(grp = group_indices(., id.x)) %>%
    gather(i, id, -grp) %>%
    distinct(grp, id) %>%
    arrange(grp, id)

# here we see that many duplicates have different alleles for some loci
# potential genotyping errors
# especially at the D1S
dups %>%
    left_join(profilesdf) %>%
    distinct(grp, marker, h1, h2) %>% 
    add_count(grp, marker) %>% 
    filter(n>1) %>%
    arrange(grp, marker) %>%
    distinct(grp, marker) %>%
    count(marker, sort = TRUE) %>% print(n=Inf)

#dups2 <- read_tsv("./sharetable_all.tsv") %>%
#    select(1:2) %>%
#    rowid_to_column("pair") %>%
#    gather(i, id, -pair) %>%
#    select(-i) %>%
#    arrange(pair, id)
#
#dups_list <- split(dups2, dups2$pair)
#
#pairslist <- map(dups_list, ~left_join(., dups2, by = "id") %>% distinct(pair.x, pair.y)) %>%
#    bind_rows() %>%
#    mutate(pa = pmin(pair.x, pair.y),
#	   pb = pmax(pair.x, pair.y)) %>%
#    distinct(pa, pb) %>%
#    filter(pa != pb)


# calculate allele frequencies
dups_to_calc <- dups %>%
    left_join(profilesdf) %>%
    distinct(grp, marker, h1, h2) %>% 
    add_count(grp, marker) %>% 
    filter(n == 1) %>%
    select(-grp, -n)

non_dups_to_calc <- profilesdf %>%
    filter(! id %in% dups$id) %>%
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

print(not_in_parents, n = Inf)


# Update allele freqs
# alleles with f > 5/2N or alleles not observed in parents will get f = 5/2N
updated_freq_df <- bind_rows(freq_df, not_in_parents) %>%
    group_by(marker) %>%
    mutate(f = ifelse(f < 5/sum(n, na.rm = TRUE) | is.na(f), 5/sum(n, na.rm = TRUE), f),
	   adj_f = f/sum(f)) %>%
    ungroup() %>%
    arrange(marker, allele)


