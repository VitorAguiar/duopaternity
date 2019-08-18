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
                         momid = c(NA,NA),
                         sex = c("female", "male"))

ped2 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c(NA, NA),
                         momid = c(NA,NA),
                         sex = c("female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../../allele_frequencies/allele_frequencies.tsv")

trios <- trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_exc <- read_tsv("../trios/trios_cpi.tsv") %>%
    filter(n_loci >= 18, n_exclusions >= 4) %>%
    select(case_no, trio) %>%
    inner_join(trios)

all_loci <- sort(unique(trios_exc$marker))
familias_loci <- map(all_loci, make_familias_locus, freqs = freqs)

duos <- trios_exc %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi),
           exclusion = as.integer(pi == 0))

duos_cpi <- duos %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(adj_pi),
              n_exclusions = sum(exclusion)) %>%
    ungroup() 

write_tsv(duos, "duos_pi_r001.tsv")
write_tsv(duos_cpi, "duos_cpi_r001.tsv")

duos_inc <- duos_cpi %>%
    filter(n_exclusions < 4 & cpi >= 10000)


