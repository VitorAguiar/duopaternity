library(tidyverse)

profilesdf <- read_tsv("../input_data/integrated_data.tsv") %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    gather(id, allele, -(case_no:marker)) %>%
    separate(id, c("ind", "h"), sep = "_") %>%
    unite(id, c("case_no", "ind"), sep = "_") %>%
    mutate(h = paste0("h", h)) %>%
    spread(h, allele)

profiles_groups <- profilesdf %>%
    group_split(marker) %>%
    map_df(~mutate(., g = group_indices(., h1, h2))) %>%
    select(id, marker, g)%>%
    arrange(marker, g, id)

write_tsv(profiles_groups, "profile_groups.tsv")
