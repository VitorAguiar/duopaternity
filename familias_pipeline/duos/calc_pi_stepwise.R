library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs, mutation) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
                  allelenames = freqs$allele[freqs$marker == locus], 
                  name = locus,
                  maleMutationModel = "Stepwise",
                  femaleMutationModel = "Equal",
                  maleMutationRate = mutation$r[mutation$marker == locus],
                  femaleMutationRate = 0,
                  maleMutationRate2 = 0.001,
                  femaleMutationRate2 = 0,
                  MutationRange = 0.1)
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

mutation_rates <- "../../input_data/aabb_male_mutationrates.tsv" %>%
    read_tsv() %>%
    filter(marker %in% freqs$marker) %>%
    mutate(r = ifelse(is.na(r), max(r, na.rm = TRUE), r))

all_loci <- sort(unique(freqs$marker))
familias_loci <- map(all_loci, make_familias_locus, freqs = freqs, mutation = mutation_rates)

str_colnames <- paste(rep(all_loci, each = 2), 1:2, sep = ".")

simul_duos <- 
    read_tsv("./duos_familias_format.tsv", 
	     col_types = cols(case_no = "i", trio = "c", person = "c", .default = "d")) %>%
    select(1:3, str_colnames)

simul_duos_long <- simul_duos %>%
    pivot_longer(-(1:3), names_to = c("marker", "h"), names_pattern = "(.+)\\.(.)", values_to = "allele") %>%
    drop_na() %>%
    unite("tocol", c("person", "h"), sep = "_") %>%
    pivot_wider(names_from = tocol, values_from = allele) %>%
    rename(ch_1 = child_1, ch_2 = child_2, af_1 = AF_1, af_2 = AF_2) %>%
    mutate(exclusion = as.integer(ch_1 != af_1 & ch_1 != af_2 & ch_2 != af_1 & ch_2 != af_2)) 

#parallelize as in trios analysis

duos_pi <- simul_duos %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
    ungroup()

duos_cpi <- duos_pi %>%
    left_join(simul_duos_long, by = c("case_no", "trio", "marker")) %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(pi),
              n_exclusions = sum(exclusion, na.rm = TRUE)) %>%
    ungroup()

duos_inclusion <- duos_cpi %>% 
    filter(n_exclusions < 4, cpi >= 10000) 

write_tsv(duos, "duos_stewise_pi.tsv")
write_tsv(duos_cpi, "duos_stepwise_cpi.tsv")
write_tsv(duos_inclusion, "duos_stepwise_falseinclusions.tsv")
