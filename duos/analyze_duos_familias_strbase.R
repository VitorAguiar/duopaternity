library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs, mutation) {

    FamiliasLocus(frequencies = freqs$f[freqs$marker == locus],
		  allelenames = freqs$allele[freqs$marker == locus], 
		  name = locus,
		  MutationModel = "Stepwise",
		  MutationRate = mutation$rate[mutation$marker == locus],
		  MutationRate2 = 0.00001,
		  MutationRange = 0.1)
}

calc_cpi <- function(df_profiles, loci, pedigrees) {
    
    datamatrix <- df_profiles %>%
	select(-case_no, -trio) %>%
	as.data.frame() %>%
	column_to_rownames("person")

    result <- FamiliasPosterior(pedigrees, loci, datamatrix, ref = 2)

    tibble(cpi = result$LR[["isFather"]])
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

freqs <- read_tsv("./allele_frequency.tsv")

mutation_rates <- read_tsv("./strbase_mutation_rates.tsv") %>%
    arrange(marker)

duos <- read_tsv("./exclusion_trios.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)


# Identifiler
identifiler_loci <- sort(readLines("./identifiler_loci.txt"))

familias_ident_loci <- map(identifiler_loci, make_familias_locus, 
			   freqs = freqs, mutation = mutation_rates)

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
   
# PP16
pp16_loci <- sort(readLines("./pp16_loci.txt"))

familias_pp16_loci <- map(pp16_loci, make_familias_locus, 
			   freqs = freqs, mutation = mutation_rates)

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


