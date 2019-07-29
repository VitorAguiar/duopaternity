library(tidyverse)

profiles <- read_tsv("../input_data/integrated_data.tsv") %>%
    filter(trio == "M1_F1_SP1") %>%
    select(-trio, -ch_1, -ch_2) %>%
    gather(hap, allele, m_1:af_2) %>%
    separate(hap, c("ind", "h"), sep = "_") %>%
    mutate(h = paste0("h", h)) %>%
    spread(h, allele) %>%
    unite(sampleid, c("case_no", "ind"), sep = "_") %>%
    mutate(id = group_indices(., sampleid)) %>%
    arrange(id, marker) %>%
    select(id, everything())

profiles_groups <- profiles %>%
    select(id, marker, h1, h2) %>%
    group_by(marker) %>%
    group_map(~mutate(.x, gid = group_indices(.x, h1, h2)), keep = TRUE) %>%
    bind_rows() %>%
    select(id, marker, gid) %>%
    arrange(marker, gid)

id_keys  <- distinct(profiles, id, sampleid)

write_tsv(profiles_groups, "profile_groups.tsv")
write_tsv(id_keys, "id_to_sampleid.tsv")
