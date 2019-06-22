library(tidyverse)

dat <- read_tsv("./integrated_data.tsv") %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2)

csf <- filter(dat, marker == "CSF1PO") %>%
    select(-marker) %>%
    gather(ind, allele, m_1:af_2) %>%
    separate(ind, c("id", "hap"), sep = "_") %>%
    arrange(case_no, id, hap)

csf_table <- crossing(csf, csf) %>%
    filter(!(case_no == case_no1 & id == id1))
