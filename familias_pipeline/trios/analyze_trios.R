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
	gather(hap, allele, 4:9) %>%
	separate(hap, c("person", "h"), sep = "_") %>%
	mutate(person = recode(person, "m" = "mother", "ch" = "child", "af" = "AF")) %>%
	unite("m", c("marker", "h"), sep = ".") %>%
	spread(m, allele)
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
	select(marker, pi = isFather)

    pi_df
}


ped1 <- FamiliasPedigree(id = c("mother", "child", "AF"), 
			 dadid = c(NA, "AF", NA),
			 momid = c(NA, "mother",NA),
			 sex = c("female", "female", "male"))

ped2 <- FamiliasPedigree(id = c("mother", "child", "AF"), 
			 dadid = c(NA, NA, NA),
			 momid = c(NA, "mother", NA),
			 sex = c("female", "female", "male"))

mypedigrees <- list(isFather = ped1, unrelated = ped2)

freqs <- read_tsv("../../input_data/allele_frequency.tsv")

all_loci <- sort(readLines("../../input_data/loci.txt"))
familias_all_loci <- map(all_loci, make_familias_locus, freqs = freqs)

codis_loci <- readLines("../../input_data/codis_loci.txt")
codisplus_loci <- readLines("../../input_data/codisplus_loci.txt")
ident_loci <- readLines("../../input_data/identifiler_loci.txt")
pp16_loci <- readLines("../../input_data/pp16_loci.txt")

trios <- read_tsv("../../input_data/integrated_data.tsv")

trios_pi <- trios %>%
    format_data() %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_all_loci, pedigrees = mypedigrees)) %>%
    ungroup()

trios_df <- left_join(trios, trios_pi, by = c("case_no", "trio", "marker")) %>%
    mutate(pi_adj = ifelse(pi == 0, 0.001, pi),
	   exclusion = as.integer(pi == 0))

write_tsv(trios_df, "./trios_pis.tsv") 

trios_cpi_original <- trios_df %>%
    group_by(case_no, trio) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(exclusion),
	      cpi = prod(pi_adj)) %>%
    ungroup() 

trios_cpi_codis <- trios_df %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(exclusion),
	      cpi = prod(pi_adj)) %>%
    ungroup()

trios_cpi_codisplus <- trios_df %>%
    filter(marker %in% codisplus_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 20) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(exclusion),
	      cpi = prod(pi_adj)) %>%
    ungroup()

trios_cpi_ident <- trios_df %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(exclusion),
	      cpi = prod(pi_adj)) %>%
    ungroup() 

trios_cpi_pp16 <- trios_df %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(exclusion),
	      cpi = prod(pi_adj)) %>%
    ungroup()

trios_sets <- 
    list(original = trios_cpi_original,
	 codis = trios_cpi_codis,
	 codisplus = trios_cpi_codisplus,
	 identifiler = trios_cpi_ident,
	 pp16 = trios_cpi_pp16) %>%
    bind_rows(.id = "marker_set") %>%
    arrange(case_no, trio, desc(n_loci), marker_set)

trios_exclusions <- trios_sets %>%
    filter(n_exclusions >= 3, cpi < 10000)


# total trios
total_trios <- 
    left_join(trios_sets %>% count(marker_set),
	      trios_exclusions %>% count(marker_set), by = c("marker_set")) %>%
    rename(total = n.x, exclusion = n.y) %>%
    mutate(marker_set = factor(marker_set, levels = c("original", "codis", "codisplus", "identifiler", "pp16"))) %>%
    arrange(marker_set)

write_tsv(total_trios, "./total_trios.tsv")

# separate sets
trios_exclusion_original <- trios_df %>%
    group_by(case_no, trio) %>%
    filter(sum(exclusion) >= 3) %>%
    ungroup() %>%
    select(-pi, -pi_adj, -exclusion)

trios_exclusion_codis <- trios_df %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18, sum(exclusion) >= 3) %>%
    ungroup() %>%
    select(-pi, -pi_adj, -exclusion)

trios_exclusion_codisplus <- trios_df %>%
    filter(marker %in% codisplus_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 20, sum(exclusion) >= 3) %>%
    ungroup() %>%
    select(-pi, -pi_adj, -exclusion)

trios_exclusion_ident <- trios_df %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(exclusion) >= 3) %>%
    ungroup() %>%
    select(-pi, -pi_adj, -exclusion)

trios_exclusion_pp16 <- trios_df %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(exclusion) >= 3) %>%
    ungroup() %>%
    select(-pi, -pi_adj, -exclusion)

write_tsv(trios_exclusion_original, "./trios_exclusion_original.tsv") 
write_tsv(trios_exclusion_codis, "./trios_exclusion_codis.tsv") 
write_tsv(trios_exclusion_codisplus, "./trios_exclusion_codisplus.tsv") 
write_tsv(trios_exclusion_ident, "./trios_exclusion_ident.tsv") 
write_tsv(trios_exclusion_pp16, "./trios_exclusion_pp16.tsv") 

