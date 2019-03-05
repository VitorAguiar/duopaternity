library(tidyverse)

abi <- read_tsv("./abi_filtered.tsv")
mega <- read_tsv("./megabace_filtered.tsv")

i <- inner_join(select(abi, 1:3), select(mega, 1:3))

mis <- anti_join(inner_join(abi, i), inner_join(mega, i)) %>%
  distinct(case_no, trio)

abi_clean <- anti_join(abi, mis)
mega_clean <- anti_join(mega, mis)

abi_sub <- filter(abi_clean, !case_no %in% mega_clean$case_no)

final_df <- bind_rows(mega_clean, abi_sub) %>%
  arrange(case_no, trio, marker)

mother_exclusions <- final_df %>%
  filter(m_1 != ch_1 & m_1 != ch_2 &
	 m_2 != ch_1 & m_2 != ch_2) %>%
  distinct(case_no, trio)

final_clean <- anti_join(final_df, mother_exclusions)

write_tsv(final_clean, "./integrated_data.tsv")
