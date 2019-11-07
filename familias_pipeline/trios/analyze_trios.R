library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
		  allelenames = freqs$allele[freqs$marker == locus], 
		  name = locus,
		  MutationModel = "Equal",
		  MutationRate = 0)
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

i_batch <- as.integer(commandArgs(TRUE))[1]

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

str_colnames <- paste(rep(all_loci, each = 2), 1:2, sep = ".")

trios <- 
    read_tsv("./trios_familias_format.tsv", 
	     col_types = cols(case_no = "i", trio = "c", "person" = "c", .default = "d")) %>%
    filter(batch == i_batch) %>%
    select(case_no, trio, person, str_colnames)

trios_pi <- trios %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_all_loci, pedigrees = mypedigrees)) %>%
    ungroup()

trios_res <- trios_pi %>%
    group_by(case_no, trio) %>%
    summarise(n_loci = n(), 
              n_exclusions = sum(pi == 0)) %>%
    ungroup() 

write_tsv(trios_res, sprintf("./trios_results_batch%d.tsv", i_batch))
