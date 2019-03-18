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
	ungroup()

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

loci <- readLines("../input_data/loci.txt")
loci <- sort(loci[loci %in% mutation_rates$marker])
familias_loci <- 
    map(loci, make_familias_locus, freqs = freqs, mutation = mutation_rates)

CHUNK <- commandArgs(TRUE)[1]

data_in <- paste0("./data/simul_chunk", CHUNK, ".tsv")
data_out <- paste0("./results/pi_strbase_chunk", CHUNK, ".tsv")

prof_df <- read_tsv(data_in) %>% 
    filter(marker %in% loci) %>%
    arrange(case_no, marker) 

exc_df <- prof_df %>%
    mutate(exclusion = as.integer(ch_1 != af_1 & ch_2 != af_2 & ch_1 != af_2 & ch_2 != af_1)) %>%
    select(case_no, marker, exclusion)

pi_df <- prof_df %>%
    format_data() %>%
    apply_calc_pi(familias_loci, mypedigrees)

out_df <- left_join(pi_df, exc_df, by = c("case_no", "marker"))

write_tsv(out_df, data_out)

