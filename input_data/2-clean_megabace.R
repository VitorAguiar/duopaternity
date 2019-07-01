library(tidyverse)

data_mega <- 
    read_lines("./old_data/dadosDNA", skip = 5) %>%
    split(cumsum(grepl("^\\d{5,}", .))) %>%
    setNames(map(., 1)) %>%
    map(~.[-1]) %>%
    map(~.[!grepl("amel", ., ignore.case = TRUE)]) %>% 
    map(~sub("^;;;;;", "", .)) %>%
    map(~.[grepl("M1\\|F\\d\\|SP1;", .)]) %>%
    keep(~length(.) >= 15L) %>%
    enframe("info", "genos") %>%
    separate(info, c("case_no", "date", "case_type", "state", "subject_info"), sep = ";") %>%
    mutate(case_no = as.integer(case_no), 
	   case_type = toupper(case_type)) %>%
    filter(case_no >= 101945L, case_type %in% c("", "TRIO"))

data_genos <- data_mega %>%
    select(case_no, genos) %>%
    unnest() %>%
    extract(genos, c("trio", "genos", "marker"), "([^;]+);(.*);([^;]+)") %>%
    mutate(marker = gsub(" ", "", marker)) %>%
    extract(genos, c("g_M", "g_F", "g_SP", "PI"), "^([0-9.]*;[0-9.]*);([0-9.]*;[0-9.]*);([0-9.]*;[0-9.]*);(.*)") %>%
    mutate(trio = gsub("\\|", "_", trio)) %>%
    separate(g_M, c("a.1_M", "a.2_M"), sep = ";") %>%
    separate(g_F, c("a.1_F", "a.2_F"), sep = ";") %>%
    separate(g_SP, c("a.1_SP", "a.2_SP"), sep = ";") %>%
    filter(a.1_M != "", a.2_M != "", a.1_F != "", a.2_F != "", 
	   a.1_SP != "", a.2_SP != "") %>% 
    mutate_at(vars(a.1_M:a.2_SP), as.numeric) %>%
    drop_na() %>%
    mutate(m_1 = pmin(a.1_M, a.2_M), m_2 = pmax(a.1_M, a.2_M),
	   ch_1= pmin(a.1_F, a.2_F), ch_2 = pmax(a.1_F, a.2_F),
	   af_1 = pmin(a.1_SP, a.2_SP), af_2 = pmax(a.1_SP, a.2_SP)) %>%
    select(case_no, trio, marker, m_1, m_2, ch_1, ch_2, af_1, af_2) %>%
    distinct()

data_genos_uniq <- data_genos %>% 
    add_count(case_no, trio, marker) %>%
    group_by(case_no, trio) %>%
    filter(all(n == 1L)) %>%
    ungroup() %>%
    select(-n) %>%
    gather(info, allele, -case_no, -trio, -marker) %>%
    group_by(case_no, trio) %>%
    filter(!all(allele %in% c(0, 7, 8, 9, 99))) %>%
    group_by(case_no, trio, marker) %>%
    filter(!any(allele > 99)) %>%
    ungroup() %>%
    spread(info, allele) %>%
    add_count(case_no, trio) %>%
    filter(n >= 15L) %>%
    select(case_no, trio, marker, m_1, m_2, ch_1, ch_2, af_1, af_2)

write_tsv(data_genos_uniq, "./megabace_filtered.tsv")
