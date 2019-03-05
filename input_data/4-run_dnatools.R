library(tidyverse)
library(DNAtools)

integrated_data <- read_tsv("./integrated_data.tsv") 

profiles_df <- integrated_data %>%
    distinct(case_no, marker, .keep_all = TRUE) %>%
    select(-trio, -ch_1, -ch_2) %>%
    gather(id, allele, m_1:af_2) %>%
    separate(id, c("subject_id", "h"), sep = "_") %>%
    unite(locus_id, c("marker", "h"), sep = ".") %>%
    unite(id, c("case_no", "subject_id"), sep = ".") %>%
    spread(locus_id, allele)

profiles_df[is.na(profiles_df)] <- 0

write_tsv(profiles_df, "./profiles.tsv")

dbres <- dbCompare(profiles_df, hit = 20, wildcard = TRUE, threads = 16)

saveRDS(dbres, "./dnatools_results.rds")
