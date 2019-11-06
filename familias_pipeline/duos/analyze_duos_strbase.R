library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs, mutation) {

    FamiliasLocus(frequencies = freqs$adj_f[freqs$marker == locus],
                  allelenames = freqs$allele[freqs$marker == locus], 
                  name = locus,
                  MutationModel = "Stepwise",
		  femaleMutationModel = "Equal",
                  MutationRate = mutation$rate[mutation$marker == locus],
		  femaleMutationRate = 0,
                  MutationRate2 = 0.001,
		  femaleMutationRate2 = 0,
                  MutationRange = 0.5,
		  femaleMutationRange = 0)
}

format_data <- function(dat) {

    dat %>%
        select(case_no, trio, marker, ch_1:af_2) %>%
        gather(hap, allele, -case_no, -trio, -marker) %>%
        separate(hap, c("person", "h"), sep = "_") %>%
        mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
        unite("m", c("marker", "h"), sep = ".") %>%
        spread(m, allele)
}

calc_pi <- function(df_profiles, loci, pedigrees) {
    
    datamatrix <- df_profiles %>%
    	select(-case_no, -trio) %>%
    	as.data.frame() %>%
    	column_to_rownames("person")

    result <- FamiliasPosterior(pedigrees, loci, datamatrix, ref = 2)

    pi_df <- result$LRperMarker[, "isFather", drop = FALSE] %>%
    	as.data.frame() %>%
    	rownames_to_column("marker") %>%
    	select(marker, pi = isFather)

    pi_df
}

ped1 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c("AF", NA),
                         momid = c(NA, NA),
                         sex = c("female", "male"))

ped2 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c(NA, NA),
                         momid = c(NA, NA),
                         sex = c("female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../../allele_frequencies/allele_frequencies.tsv")

mutation_rates <- "../../input_data/aabb_mutation_rates.tsv" %>%
    read_tsv() %>%
    filter(marker %in% freqs$marker) %>%
    complete(marker = freqs$marker, fill = list(rate = 0))

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_exc <- read_tsv("../trios/trios_cpi.tsv") %>%
    filter(n_loci >= 18, n_exclusions >= 4) %>%
    select(case_no, trio) %>%
    inner_join(trios)

all_loci <- sort(unique(trios_exc$marker))
familias_loci <- map(all_loci, make_familias_locus, freqs = freqs, mutation = mutation_rates)

duos_exc <- trios_exc %>%
    mutate(exclusion = as.integer(ch_1 != af_1 & ch_1 != af_2 & ch_2 != af_1 & ch_2 != af_2)) %>%
    select(case_no, trio, marker, exclusion)

duos <- trios_exc %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    inner_join(duos_exc, ., by = c("case_no", "trio", "marker")) %>%
    mutate(adj_pi = ifelse(exclusion == 1 & pi == 0, 0.001, pi)) %>%
    select(case_no, trio, marker, pi, adj_pi, exclusion)
    
duos_cpi <- duos %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(adj_pi),
              n_exclusions = sum(exclusion)) %>%
    ungroup()

duos_inclusion <- duos_cpi %>% 
    filter(n_exclusions < 4, cpi >= 10000) 

write_tsv(duos, "duos_pi_stewise.tsv")
write_tsv(duos_cpi, "duos_cpi_stepwise.tsv")
write_tsv(duos_inclusion, "duos_inclusion_stepwise.tsv")
