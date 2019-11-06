library(Familias)
library(tidyverse)
library(parallel)

make_familias_locus <- function(locus, freqs) {
    
    FamiliasLocus(frequencies = freqs$adj_f[freqs$marker == locus],
                  allelenames = freqs$allele[freqs$marker == locus], 
                  name = locus,
                  MutationModel = "Equal",
                  MutationRate = 0)
}

format_data <- function(dat) {
    
    dat %>%
        gather(hap, allele, -case_no, -marker) %>%
        separate(hap, c("person", "h"), sep = "_") %>%
        mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
        unite("m", c("marker", "h"), sep = ".") %>%
        spread(m, allele)
} 

calc_pi <- function(df_profiles, loci, pedigrees) {
    
    datamatrix <- df_profiles %>%
        select(-case_no) %>%
        as.data.frame() %>%
        column_to_rownames("person")
    
    result <- FamiliasPosterior(pedigrees, loci, datamatrix, ref = 2)
    
    pi_df <- result$LRperMarker[, "isFather", drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("marker") %>%
        select(marker, pi = isFather)
    
    pi_df
}

apply_calc_pi <- function(dat) {
    dat %>%
        format_data() %>%
        group_by(case_no) %>%
        do(calc_pi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
        ungroup() %>%
        inner_join(duos_exc, ., by = c("case_no", "marker")) %>%
        mutate(adj_pi = ifelse(exclusion == 1 & pi == 0, 0.001, pi)) %>%
        select(case_no, marker, pi, adj_pi, exclusion)
}


ped1 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c("AF", NA),
                         momid = c(NA,NA),
                         sex = c("female", "male"))

ped2 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c(NA, NA),
                         momid = c(NA,NA),
                         sex = c("female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../allele_frequencies/allele_frequencies.tsv")

duos <- data.table::fread("./data/simulated_duos.tsv") %>%
    as_tibble()

all_loci <- sort(unique(duos$marker))
familias_loci <- map(all_loci, make_familias_locus, freqs = freqs)

duos_exc <- duos %>%
    mutate(exclusion = as.integer(ch_1 != af_1 & ch_1 != af_2 & ch_2 != af_1 & ch_2 != af_2)) %>%
    select(case_no, marker, exclusion)

duos_list <- duos %>%
    mutate(i = ntile(case_no, 10)) %>%
    split(.$i) %>%
    map(~select(., -i))

duos_pi <- mclapply(duos_list, apply_calc_pi, mc.cores = 10) %>%
    bind_rows()

duos_cpi <- duos_pi %>%
    group_by(case_no) %>%
    summarise(cpi = prod(adj_pi),
              n_exclusions = sum(exclusion)) %>%
    ungroup() 

duos_inclusion <- duos_cpi %>%
    filter(n_exclusions < 4 & cpi >= 10000)

write_tsv(duos_pi, "./duos_pi_r001.tsv")
write_tsv(duos_cpi, "./duos_cpi_r001.tsv")
write_tsv(duos_inclusion, "./duos_inclusion_r001.tsv")
