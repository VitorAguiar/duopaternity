library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
		  allelenames = freqs$allele[freqs$marker == locus], 
		  name = locus,
		  MutationModel = "Equal",
		  MutationRate = 0)
}

format_data <- function(duos_data) {

    duos_data %>%
	gather(hap, allele, 3:6) %>%
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

apply_calc_pi <- function(duos_data, familias_loci, familias_pedigrees) {
    
    duos_data %>%
	group_by(case_no) %>%
	do(calc_pi(., loci = familias_loci, pedigrees = familias_pedigrees)) %>%
	ungroup() %>%
	mutate(adj_pi = ifelse(pi == 0, 0.001, pi))
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

freqs <- read_tsv("../input_data/allele_frequency.tsv")

loci <- sort(readLines("../input_data/loci.txt"))
familias_loci <- map(loci, make_familias_locus, freqs = freqs)


CHUNK <- commandArgs(TRUE)[1]

data_in <- paste0("./data/simul_chunk", CHUNK, ".tsv")
data_out <- paste0("./results/pi_r001_chunk", CHUNK, ".tsv") 

pi_df <- read_tsv(data_in) %>%
    format_data() %>%
    apply_calc_pi(familias_loci, mypedigrees)

write_tsv(pi_df, data_out)

