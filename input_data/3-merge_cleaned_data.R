library(tidyverse)

strinfo <- read_tsv("./str_info.tsv")

abi <- read_tsv("./abi_filtered.tsv")
mega <- read_tsv("./megabace_filtered.tsv")

abilong <- abi %>% 
    gather(hap, allele, m_1:af_2)

megalong <- mega %>% 
    gather(hap, allele, m_1:af_2)

notcompatible <- inner_join(abilong, megalong, by = c("case_no", "trio", "marker", "hap")) %>%
    filter(allele.x != allele.y) %>%
    distinct(case_no, trio, marker)

abi_clean <- anti_join(abi, notcompatible)
mega_clean <- anti_join(mega, notcompatible)

merge_df <- bind_rows(mega_clean, abi_clean) %>%
    distinct(case_no, trio, marker, .keep_all = TRUE) %>%
    arrange(case_no, trio, marker) %>%
    add_count(case_no, trio) %>%
    filter(n >= 15L) %>%
    select(-n)

mother_exclusions <- merge_df %>%
    filter(m_1 != ch_1 & m_1 != ch_2 & m_2 != ch_1 & m_2 != ch_2) %>%
    distinct(case_no, trio)

merge_clean <- anti_join(merge_df, mother_exclusions)

final_loci <- merge_clean %>% 
    count(marker) %>%
    filter(n >= 2000) %>%
    pull(marker)

final_merge_data <- merge_clean %>%
    filter(marker %in% final_loci) %>%
    add_count(case_no, trio) %>%
    filter(n >= 15L) %>%
    select(-n)

final_clean <- final_merge_data %>% 
    gather(hap, allele, m_1:af_2) %>%
    mutate(microvariant = allele %% 1L * 10L) %>%
    left_join(strinfo, by = c("marker" = "locus")) %>%
    filter(microvariant < repeats) %>%
    select(case_no:allele) %>%
    spread(hap, allele) %>%
    drop_na() %>%
    select(case_no, trio, marker, m_1, m_2, ch_1, ch_2, af_1, af_2) %>%
    add_count(case_no, trio) %>%
    filter(n >= 15L) %>%
    select(-n)

write_tsv(final_clean, "./integrated_data.tsv")
