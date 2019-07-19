library(tidyverse)

abi <- read_tsv("./abi_filtered.tsv")
mega <- read_tsv("./megabace_filtered.tsv")

i <- inner_join(select(abi, 1:3), select(mega, 1:3))

mis <- anti_join(inner_join(abi, i), inner_join(mega, i)) %>%
    distinct(case_no, trio, marker)

abi_clean <- anti_join(abi, mis)
mega_clean <- anti_join(mega, mis)

abi_sub <- filter(abi_clean, !case_no %in% mega_clean$case_no)

merge_df <- bind_rows(mega_clean, abi_sub) %>%
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

write_tsv(final_merge_data, "./integrated_data.tsv")
