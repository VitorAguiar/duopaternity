library(Familias)
library(tidyverse)

make_familias_locus <- function(locus, freqs, mutation) {

    FamiliasLocus(frequencies = freqs$adj_f[freqs$marker == locus],
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
                         momid = c(NA, NA),
                         sex = c("female", "male"))

ped2 <- FamiliasPedigree(id = c("child", "AF"), 
                         dadid = c(NA, NA),
                         momid = c(NA, NA),
                         sex = c("female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../../allele_frequencies/allele_frequencies.tsv")

mutation_rates <- "../../input_data/aabb_mutation_rates.tsv" %>%
    read_tsv()

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_exc <- read_tsv("../trios/trios_cpi.tsv") %>%
    filter(n_loci >= 18, n_exclusions >= 4) %>%
    select(case_no, trio) %>%
    inner_join(trios)

all_loci <- sort(unique(trios_exc$marker))


all_loci[! all_loci %in% mutation_rates$marker]


familias_ident_loci <- map(all_loci, make_familias_locus, freqs = freqs,
                           mutation = mutation_rates)

duos_ident <- read_tsv("../trios/trios_exclusion_ident.tsv")

duos_ident_exc <- duos_ident %>%
    mutate(exclusion = as.integer(ch_1 != af_1 & ch_1 != af_2 & ch_2 != af_1 & ch_2 != af_2)) %>%
    select(case_no, trio, marker, exclusion)

duos_ident_pi <- duos_ident %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_pi(., loci = familias_ident_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    left_join(duos_ident_exc, by = c("case_no", "trio", "marker"))

duos_ident_cpi <- duos_ident_pi %>%
    group_by(case_no, trio) %>%
    summarise(cpi = prod(pi),
	      n_exclusions = sum(exclusion)) %>%
    ungroup()

write_tsv(duos_ident_pi, "duos_ident_pi_strbase.tsv")
write_tsv(duos_ident_cpi, "duos_ident_cpi_strbase.tsv")

duos_ident_inc <- duos_ident_cpi %>%
    filter(n_exclusions < 3 & cpi >= 10000)


