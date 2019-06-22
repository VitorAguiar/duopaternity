library(tidyverse)

opts <- commandArgs(TRUE)
profiles_df <- opts[1]
chunk       <- as.integer(opts[2])
nchunks     <- as.integer(opts[3])
outPrefix   <- opts[4]

out <- sprintf("%s_%d.tsv", outPrefix, chunk)

dat <- profiles_df %>% 
    read_tsv() %>%
    filter(trio == "M1_F1_SP1") %>%
    select(case_no, marker, m_1, m_2, af_1, af_2)

dat_long <- dat %>%
    mutate(case_no = as.character(case_no)) %>%
    gather(ind, allele, m_1:af_2) %>%
    separate(ind, c("parent", "hap"), sep = "_") %>%
    unite(id, c("case_no", "parent"), sep = "_")

subjects <- unique(dat_long$id)

chunk_size <- ceiling(length(subjects)/nchunks)
chunk_start <- chunk_size * chunk - chunk_size + 1L
chunk_end <- chunk_start + chunk_size - 1L

dat_sub <- dat_long %>%
    filter(id %in% subjects[chunk_start:chunk_end])

dat_table <- left_join(dat_sub, dat_sub, by = "marker") %>%
    drop_na() %>%
    filter(id.x != id.y)

share_table <- dat_table %>%
    mutate(i1 = as.integer(hap.x == hap.y & allele.x == allele.y),
	   i2 = as.integer(hap.x != hap.y & allele.x == allele.y)) %>%
    group_by(id.x, id.y, marker) %>%
    summarise(i1 = sum(i1), i2 = sum(i2)) %>%
    ungroup() %>%
    mutate(shr = pmax(i1, i2)) %>%
    select(-i1, -i2) %>%
    group_by(id.x, id.y) %>%
    summarise(total = mean(shr == 2),
	      partial = mean(shr == 1))

write_tsv(share_table, out)
