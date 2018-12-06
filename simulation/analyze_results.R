library(tidyverse)

simul_codis_001 <- read_tsv("./results/duos_codis_pis.tsv") %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi)) %>%
    ungroup()

pval_codis_001 <- simul_codis_001 %>%
    summarise(p = mean(cpi >= 10000))

simul_ident_001 <- read_tsv("./results/duos_identifiler_pis.tsv") %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi)) %>%
    ungroup()

pval_ident_001 <- simul_ident_001 %>%
    summarise(p = mean(cpi >= 10000))

simul_pp16_001 <- read_tsv("./results/duos_pp16_pis.tsv") %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi)) %>%
    ungroup()

pval_pp16_001 <- simul_pp16_001 %>%
    summarise(p = mean(cpi >= 10000))

simul_ident_base <- read_tsv("./results/duos_identifiler_pis_strbase.tsv") %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi)) %>%
    ungroup()

pval_ident_base <- simul_ident_base %>%
    summarise(p = mean(cpi >= 10000))

simul_pp16_base <- read_tsv("./results/duos_pp16_pis_strbase.tsv") %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi)) %>%
    ungroup()

pval_pp16_base <- simul_pp16_base %>%
    summarise(p = mean(cpi >= 10000))

pval_df <- bind_rows(codis_r001 = pval_codis_001,
                     identifiler_r001 = pval_ident_001,
                     pp16_r001 = pval_pp16_001,
                     identifiler_STRbase = pval_ident_base,
                     pp16_STRbase = pval_pp16_base, 
                     .id = "id") %>%
    separate(id, c("marker_set", "mutation_rate"), sep = "_")


total_cases <- read_tsv("../familias_pipeline/trios/total_trios.tsv") %>%
    select(marker_set, N = exclusion)
    
obs_inclusions_r001 <- 
    read_tsv("../familias_pipeline/duos/false_inclusions_r001.tsv") %>%
    count(marker_set) %>%
    mutate(mutation_rate = "r001") %>%
    left_join(total_cases) %>%
    mutate(obs = n/N) %>%
    select(marker_set, mutation_rate, obs)

obs_inclusions_base <- 
    read_tsv("../familias_pipeline/duos/false_inclusions_strbase.tsv") %>%
    count(marker_set) %>%
    mutate(mutation_rate = "STRbase") %>%
    left_join(total_cases) %>%
    mutate(obs = n/N) %>%
    select(marker_set, mutation_rate, obs)

obs_df <- bind_rows(obs_inclusions_r001, obs_inclusions_base)

out <- left_join(obs_df, pval_df)

write_tsv(out, "final_results.tsv")

