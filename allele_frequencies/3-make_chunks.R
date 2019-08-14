library(tidyverse)

chunk_df <- read_tsv("../input_data/integrated_data.tsv") %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2) %>%
    gather(id, allele, -(case_no:marker)) %>%
    separate(id, c("ind", "h"), sep = "_") %>%
    unite(id, c("case_no", "ind"), sep = "_") %>%
    distinct(id) %>%
    arrange(id) %>%
    mutate(ck = ntile(id, 100))

write_tsv(chunk_df, "./chunk_df.tsv")
