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


# CODIS
codis_loci <- sort(readLines("../input_data/codis_loci.txt"))
familias_codis_loci <- map(codis_loci, make_familias_locus, freqs = freqs)

duos_codis <- read_tsv("./simulated_duos.tsv") %>%
    filter(marker %in% codis_loci) %>%
    gather(hap, allele, 3:6) %>%
    separate(hap, c("person", "h"), sep = "_") %>%
    mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    spread(m, allele)

duos_codis_pis <- duos_codis %>%
    group_by(case_no) %>%
    do(calc_pi(., loci = familias_codis_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi))

write_tsv(duos_codis_pis, "./results/duos_codis_pis.tsv")


# Identifiler
identifiler_loci <- sort(readLines("../input_data/identifiler_loci.txt"))
familias_ident_loci <- map(identifiler_loci, make_familias_locus, freqs = freqs)

duos_ident <- read_tsv("./data/simulated_duos.tsv") %>%
    filter(marker %in% identifiler_loci) %>%
    gather(hap, allele, 3:6) %>%
    separate(hap, c("person", "h"), sep = "_") %>%
    mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    spread(m, allele)

duos_ident_pis <- duos_ident %>%
    group_by(case_no) %>%
    do(calc_pi(., loci = familias_ident_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi))

write_tsv(duos_ident_pis, "./results/duos_identifiler_pis.tsv")


# PP16
pp16_loci <- sort(readLines("../input_data/pp16_loci.txt"))
familias_pp16_loci <- map(pp16_loci, make_familias_locus, freqs = freqs)

duos_pp16 <- read_tsv("./data/simulated_duos.tsv") %>%
    filter(marker %in% pp16_loci) %>%
    gather(hap, allele, 3:6) %>%
    separate(hap, c("person", "h"), sep = "_") %>%
    mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    spread(m, allele)

duos_pp16_pis <- duos_pp16 %>%
    group_by(case_no) %>%
    do(calc_pi(., loci = familias_pp16_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi))

write_tsv(duos_pp16_pis, "./results/duos_pp16_pis.tsv")
