library(tidyverse)

dat <- "../input_data/integrated_data.tsv" %>% 
    read_tsv() %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2)

dat_long <- dat %>%
    mutate(case_no = format(case_no, scientific = FALSE)) %>%
    gather(ind, allele, m_1:af_2) %>%
    separate(ind, c("parent", "hap"), sep = "_") %>%
    unite(id, c("case_no", "parent"), sep = "_")

write_tsv(dat_long, "./str_parents.tsv")
