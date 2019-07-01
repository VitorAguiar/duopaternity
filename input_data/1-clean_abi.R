library(tidyverse)
  
read_data <- function(path) {

    read_csv(path, guess_max = 10000) %>%
    select(matches("sample.*name|run.*name|marker|allele", ignore.case = TRUE)) %>%
    setNames(names(.) %>% tolower() %>% trimws() %>% gsub(" ", "_", .)) %>%
    separate(sample_name, c("case_no", "subject_id"), sep = "-", extra = "drop") %>%
    mutate_all(trimws) %>%
    mutate(case_no = as.integer(case_no)) %>%
    mutate_at(vars(allele_1, allele_2), as.numeric) %>%
    drop_na() %>%
    filter(case_no >= 206807L, grepl("F\\d|M1|SP1", subject_id)) %>%
    group_by(case_no) %>%
    filter(all(c("F1", "M1", "SP1") %in% subject_id)) %>%
    ungroup() %>%
    mutate(a1 = pmin(allele_1, allele_2), a2 = pmax(allele_1, allele_2)) %>%
    select(case_no, subject_id, run_name, marker, allele_1 = a1, allele_2 = a2)
}

#loci <- readLines("./loci.txt")

kit1 <- list.files("./ngm", full.names = TRUE, pattern = "^\\d+.*csv$") %>%
    map_df(read_data)

kit2 <- list.files("./profiler", full.names = TRUE, pattern = "^\\d+.*csv$") %>%
    map_df(read_data)

mistmatch_kits <- 
    inner_join(kit1, kit2, by = c("case_no", "subject_id", "marker")) %>%
    filter(allele_1.x != allele_1.y | allele_2.x != allele_2.y) %>%
    distinct(case_no) %>%
    pull(case_no)

kit1_clean <- filter(kit1, !case_no %in% mistmatch_kits)
kit2_clean <- filter(kit2, !case_no %in% mistmatch_kits)

kits_merge <- bind_rows(kit1_clean, kit2_clean) %>%
    #filter(marker %in% loci, allele_1 < 99, allele_2 < 99) %>%
    filter(allele_1 < 99, allele_2 < 99) %>%
    drop_na() %>%
    distinct(case_no, subject_id, marker, allele_1, allele_2) %>%
    add_count(case_no, subject_id, marker) %>%
    filter(n == 1L) %>%
    select(-n)

kits_parents <- kits_merge %>%
    filter(subject_id %in% c("M1", "SP1")) %>%
    unite(genot, c("allele_1", "allele_2"), sep = "_") %>%
    spread(subject_id, genot) %>%
    separate(M1, c("m_1", "m_2"), sep = "_") %>%
    separate(SP1, c("af_1", "af_2"), sep = "_")

kits_children <- kits_merge %>%
    filter(grepl("F\\d", subject_id)) %>%
    mutate(trio = paste0("M1_", subject_id, "_SP1")) %>%
    rename(ch_1 = allele_1, ch_2 = allele_2)

kits_out <- 
    left_join(kits_parents, kits_children, by = c("case_no", "marker")) %>%
    select(case_no, trio, marker, m_1, m_2, ch_1, ch_2, af_1, af_2) %>%
    drop_na() %>%
    add_count(case_no, trio) %>%
    filter(n >= 15L) %>%
    select(-n)

write_tsv(kits_out, "./abi_filtered.tsv")
