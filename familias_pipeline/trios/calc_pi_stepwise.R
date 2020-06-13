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
                  maleMutationRate2 = 0.000001,
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

ped1 <- FamiliasPedigree(id = c("mother", "child", "AF"), 
                         dadid = c(NA, "AF", NA),
                         momid = c(NA, "mother",NA),
                         sex = c("female", "female", "male"))

ped2 <- FamiliasPedigree(id = c("mother", "child", "AF"), 
                         dadid = c(NA, NA, NA),
                         momid = c(NA, "mother", NA),
                         sex = c("female", "female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../../allele_frequencies/allele_frequencies.tsv")

mutation_rates <- "../../input_data/aabb2008_paternal_mutationrates.tsv" %>%
    read_tsv() %>%
    filter(marker %in% freqs$marker) %>%
    mutate(r = ifelse(is.na(r), max(r, na.rm = TRUE), r))

all_loci <- sort(unique(freqs$marker))
familias_all_loci <- map(all_loci, make_familias_locus, freqs = freqs, mutation = mutation_rates)

str_colnames <- paste(rep(all_loci, each = 2), 1:2, sep = ".")

trios <- read_tsv("../../input_data/integrated_data.tsv") %>%
    filter(marker %in% freqs$marker) %>%
    group_by(case_no, trio) %>%
    filter(n() > 18) %>%
    ungroup()

trios_familias_format <- trios %>%
    pivot_longer(m_1:af_2, 
                 names_to = c("person", "h"), 
                 names_pattern = "(.+)_(.)", 
                 values_to = "allele") %>%
    mutate(person = recode(person, "m" = "mother", "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    pivot_wider(names_from = m, values_from = allele) %>%
    select(case_no, trio, person, str_colnames)

trios_pi <- trios_familias_format %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_all_loci, pedigrees = mypedigrees)) %>%
    ungroup()

trios_res <- trios_pi %>%
    group_by(case_no, trio) %>%
    summarise(n_loci = n(),
              cpi = prod(pi)) %>%
    ungroup() 

write_tsv(trios_pi, "./trios_pi_sw.tsv")
write_tsv(trios_res, "./trios_cpi_sw.tsv")
