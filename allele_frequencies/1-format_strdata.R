library(tidyverse)

str_parents <- "../input_data/integrated_data.tsv" %>% 
    read_tsv() %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    mutate(case_no = format(case_no, scientific = FALSE)) %>%
    pivot_longer(m_1:af_2, 
                 names_to = c("parent", "hap"),
                 names_pattern = "(.+)_(.+)",
                 values_to = "allele") %>%
    unite(id, c("case_no", "parent"), sep = "_") 

profile_groups <- str_parents %>%
    mutate(hap = paste0("h", hap)) %>%
    pivot_wider(names_from = hap, values_from = allele) %>%
    mutate(g = group_indices(., marker, h1, h2)) %>%
    select(id, marker, g) %>%
    arrange(marker, g, id)

chunk_df <- str_parents %>%
    distinct(id) %>%
    arrange(id) %>%
    mutate(ck = ntile(id, 100))

out <- left_join(profile_groups, chunk_df) %>%
    arrange(ck, id, marker)

write_tsv(out, "profile_groups.tsv")
