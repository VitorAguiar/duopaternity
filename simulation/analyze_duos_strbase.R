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

mutation_rates <- read_tsv("../input_data/strbase_mutation_rates.tsv") %>%
    arrange(marker)


kit <- commandArgs(TRUE)[1] #identifiler or pp16

loci <- sort(readLines(paste0("../input_data/", kit, "_loci.txt")))

familias_loci <- 
    map(loci, make_familias_locus, freqs = freqs, mutation = mutation_rates)

duos <- read_tsv("./data/simulated_duos.tsv") %>%
    filter(marker %in% loci) %>%
    gather(hap, allele, 3:6) %>%
    separate(hap, c("person", "h"), sep = "_") %>%
    mutate(person = recode(person, "ch" = "child", "af" = "AF")) %>%
    unite("m", c("marker", "h"), sep = ".") %>%
    spread(m, allele)

duos_pis <- duos %>%
    group_by(case_no) %>%
    do(calc_pi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(adj_pi = ifelse(pi == 0, 0.001, pi))

write_tsv(duos_pis, paste0("./results/duos_", kit, "_pis_strbase.tsv"))
