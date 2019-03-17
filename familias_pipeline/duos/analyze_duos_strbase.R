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

format_data <- function(dat) {

    dat %>%
        select(case_no, trio, marker, ch_1:af_2) %>%
        gather(hap, allele, -case_no, -trio, -marker) %>%
        separate(hap, c("person", "h"), sep = "_") %>%
        mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
        unite("m", c("marker", "h"), sep = ".") %>%
        spread(m, allele)
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

freqs <- read_tsv("../../input_data/allele_frequency.tsv")

mutation_rates <- read_tsv("../../input_data/strbase_mutation_rates.tsv") %>%
    arrange(marker)


# Identifiler
ident_loci <- sort(readLines("../../input_data/identifiler_loci.txt"))
familias_ident_loci <- map(ident_loci, make_familias_locus, freqs = freqs,
			   mutation = mutation_rates)

duos_ident <- read_tsv("../trios/trios_exclusion_ident.tsv")

duos_ident_exc <- duos_ident %>%
    group_by(case_no, trio) %>%
    summarise(n_exclusions = sum(exclusion)) %>%
    ungroup()

duos_ident_pis <- duos_ident %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_ident_loci, pedigrees = mypedigrees)) %>%
    ungroup()

duos_ident_inc <- duos_ident_pis %>%
    left_join(duos_ident_exc, by = c("case_no", "trio")) %>%
    summarise(n = sum(n_exclusions < 3 & cpi >= 10000))


# PP16
pp16_loci <- sort(readLines("../../input_data/pp16_loci.txt"))
familias_pp16_loci <- map(pp16_loci, make_familias_locus, freqs = freqs,
			   mutation = mutation_rates)

duos_pp16 <- read_tsv("../trios/trios_exclusion_pp16.tsv")

duos_pp16_exc <- duos_pp16 %>%
    group_by(case_no, trio) %>%
    summarise(n_exclusions = sum(exclusion)) %>%
    ungroup()

duos_pp16_pis <- duos_pp16 %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_pp16_loci, pedigrees = mypedigrees)) %>%
    ungroup()

duos_pp16_inc <- duos_pp16_pis %>%
    left_join(duos_pp16_exc, by = c("case_no", "trio")) %>%
    summarise(n = sum(n_exclusions < 3 & cpi >= 10000))

inclusion_df <- bind_rows(identifiler = duos_ident_inc,
                          pp16 = duos_pp16_inc, .id = "marker_set")

write_tsv(inclusion_df, "./false_inclusions_strbase.tsv")

