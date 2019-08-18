library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs) {

    FamiliasLocus(frequencies = freqs$adj_f[freqs$marker == locus],
		  allelenames = freqs$allele[freqs$marker == locus], 
		  name = locus,
		  MutationModel = "Equal",
		  MutationRate = 0)
}

format_data <- function(dat) {

    dat %>%
    	gather(hap, allele, 4:9) %>%
    	separate(hap, c("person", "h"), sep = "_") %>%
    	mutate(person = recode(person, "m" = "mother", "ch" = "child", "af" = "AF")) %>%
    	unite("m", c("marker", "h"), sep = ".") %>%
    	spread(m, allele)
}

calc_cpi <- function(df_profiles, loci, pedigrees) {
    
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

all_loci <- sort(unique(freqs$marker))
familias_all_loci <- map(all_loci, make_familias_locus, freqs = freqs)

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_pi <- trios %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_all_loci, pedigrees = mypedigrees)) %>%
    ungroup()

trios_df <- left_join(trios, trios_pi, by = c("case_no", "trio", "marker")) %>%
    mutate(pi_adj = ifelse(pi == 0, 0.001, pi),
           exclusion = as.integer(pi == 0))

trios_cpi <- trios_df %>%
    group_by(case_no, trio) %>%
    summarise(n_loci = n(), 
              n_exclusions = sum(exclusion),
              cpi = prod(pi_adj)) %>%
    ungroup() 

write_tsv(trios_cpi, "./trios_cpi.tsv")