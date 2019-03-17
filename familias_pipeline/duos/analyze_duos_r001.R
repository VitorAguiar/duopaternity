library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
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

freqs <- read_tsv("../../input_data/allele_frequency.tsv")


# CODIS
codis_loci <- sort(readLines("../../input_data/codis_loci.txt"))
familias_codis_loci <- map(codis_loci, make_familias_locus, freqs = freqs)

duos_codis <- read_tsv("../trios/trios_exclusion_codis.tsv") %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_codis_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi),
	   exclusion = as.integer(pi == 0))

duos_codis_cpi <- duos_codis %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup() 

write_tsv(duos_codis_cpi, "duos_codis_cpi.tsv")

duos_codis_inc <- duos_codis_cpi %>%
    filter(n_exclusions < 3 & cpi >= 10000)
    

# Identifiler
identifiler_loci <- sort(readLines("../../input_data/identifiler_loci.txt"))
familias_ident_loci <- map(identifiler_loci, make_familias_locus, freqs = freqs)

duos_ident <- read_tsv("../trios/trios_exclusion_ident.tsv") %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_ident_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi),
	   exclusion = as.integer(pi == 0))

duos_ident_cpi <- duos_ident %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(duos_ident_cpi, "duos_ident_cpi.tsv")

duos_ident_inc <- duos_ident_cpi %>%
    filter(n_exclusions < 3 & cpi >= 10000)


# PP16
pp16_loci <- sort(readLines("../../input_data/pp16_loci.txt"))
familias_pp16_loci <- map(pp16_loci, make_familias_locus, freqs = freqs)

duos_pp16 <- read_tsv("../trios/trios_exclusion_pp16.tsv") %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_pp16_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi),
	   exclusion = as.integer(pi == 0))

duos_pp16_cpi <- duos_pp16 %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(adj_pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup() 

write_tsv(duos_pp16_cpi, "duos_pp16_cpi.tsv")

duos_pp16_inc <- duos_pp16_cpi %>%
    filter(n_exclusions < 3 & cpi >= 10000)


inclusion_df <- bind_rows(codis = duos_codis_inc,
			  identifiler = duos_ident_inc,
			  pp16 = duos_pp16_inc, .id = "marker_set")

write_tsv(inclusion_df, "./false_inclusions_r001.tsv")
