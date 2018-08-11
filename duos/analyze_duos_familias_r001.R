library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
		  allelenames = freqs$allele[freqs$marker == locus], 
		  name = locus,
		  MutationModel = "Equal",
		  MutationRate = 0)
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
	mutate(pi = ifelse(isFather == 0, 0.001, isFather)) %>%
	select(marker, pi)

    pi_df %>% summarise(cpi = prod(pi))
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

duos <- read_tsv("../trios/exclusion_trios.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)


# CODIS
codis_loci <- sort(readLines("../input_data/codis_loci.txt"))
familias_codis_loci <- map(codis_loci, make_familias_locus, freqs = freqs)

duos_codis <- duos %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18) %>%
    ungroup() %>%
    arrange(case_no, trio, marker) %>%
    gather(hap, allele, 4:7) %>%
    separate(hap, c("h", "person"), sep = "_") %>%
    mutate(person = recode(person, "F" = "child", "SP" = "AF"),
	   h = sub("^a", ".", h)) %>%
    unite("m", marker:h, sep = "") %>%
    spread(m, allele)

inclusion_duos_codis <- duos_codis %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_codis_loci, pedigrees = mypedigrees)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

perc_codis <- 
    tibble(marker_set = "codis",
	   perc = (nrow(inclusion_duos_codis) / nrow(count(duos_codis, case_no, trio))) * 100)


# Identifiler
identifiler_loci <- sort(readLines("..input_data/identifiler_loci.txt"))
familias_ident_loci <- map(identifiler_loci, make_familias_locus, freqs = freqs)

duos_identifiler <- duos %>%
    filter(marker %in% identifiler_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    ungroup() %>%
    arrange(case_no, trio, marker) %>%
    gather(hap, allele, 4:7) %>%
    separate(hap, c("h", "person"), sep = "_") %>%
    mutate(person = recode(person, "F" = "child", "SP" = "AF"),
	   h = sub("^a", ".", h)) %>%
    unite("m", marker:h, sep = "") %>%
    spread(m, allele)

inclusion_duos_identifiler <- duos_identifiler %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_ident_loci, pedigrees = mypedigrees)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

perc_identifiler <-
    tibble(marker_set = "identifiler",
	   perc = (nrow(inclusion_duos_identifiler) / nrow(count(duos_identifiler, case_no, trio))) * 100)


# PP16
pp16_loci <- sort(readLines("../input_data/pp16_loci.txt"))
familias_pp16_loci <- map(pp16_loci, make_familias_locus, freqs = freqs)

duos_pp16 <- duos %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    ungroup() %>%
    arrange(case_no, trio, marker) %>%
    gather(hap, allele, 4:7) %>%
    separate(hap, c("h", "person"), sep = "_") %>%
    mutate(person = recode(person, "F" = "child", "SP" = "AF"),
	   h = sub("^a", ".", h)) %>%
    unite("m", marker:h, sep = "") %>%
    spread(m, allele)

inclusion_duos_pp16 <- duos_pp16 %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_pp16_loci, pedigrees = mypedigrees)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

perc_pp16 <- 
    tibble(marker_set = "pp16",
	   perc = (nrow(inclusion_duos_pp16) / nrow(count(duos_pp16, case_no, trio))) * 100)

inclusion_df <- list(codis = inclusion_duos_codis, 
		     identifiler = inclusion_duos_identifiler,
		     pp16 = inclusion_duos_pp16) %>%
    bind_rows(.id = "marker_set")

perc_df <- bind_rows(perc_codis, perc_identifiler, perc_pp16) 

