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

freqs <- read_tsv("../input_data/allele_frequency.tsv")

loci <- sort(readLines("../input_data/loci.txt"))
familias_loci <- map(loci, make_familias_locus, freqs = freqs)

trios <- read_tsv("../input_data/integrated_data.tsv") %>%
    setNames(sub("allele\\.", "a", names(.))) 

trios_w <- trios %>%
    gather(hap, allele, 4:9) %>%
    separate(hap, c("h", "person"), sep = "_") %>%
    mutate(person = recode(person, "M" = "mother", "F" = "child", "SP" = "AF"),
	   h = sub("^a", ".", h)) %>%
    unite("m", marker:h, sep = "") %>%
    spread(m, allele)

trios_pi <- trios_w %>%
    group_by(case_no, trio) %>%
    do(calc_cpi(., loci = familias_loci, pedigrees = mypedigrees)) %>%
    ungroup() %>%
    mutate(pi_adj = ifelse(pi == 0, 0.001, pi))

trios_df <- left_join(trios, trios_pi, by = c("case_no", "trio", "marker")) 

trios_cpi <- trios_df %>%
    group_by(case_no, trio) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(pi == 0),
	      cpi = prod(pi_adj)) %>%
    ungroup() %>%
    mutate(marker_set = "original")


perc_exclusion_per_locus <- trios_df %>%
    group_by(marker) %>%
    summarise(n = n(),
	      ne = sum(pi == 0)) %>%
    ungroup() %>%
    mutate(perc = ne/n * 100) %>%
    select(marker, perc) %>%
    arrange(desc(perc))

write_tsv(perc_exclusion_per_locus, "./perc_exclusions_per_locus.tsv")

codis_loci <- readLines("../input_data/codis_loci.txt")

trios_codis_cpi <- trios_df %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(pi == 0),
	      cpi = prod(pi_adj)) %>%
    ungroup() %>%
    mutate(marker_set = "codis")

ident_loci <- readLines("../input_data/identifiler_loci.txt")

trios_ident_cpi <- trios_df %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(pi == 0),
	      cpi = prod(pi_adj)) %>%
    ungroup() %>%
    mutate(marker_set = "identifiler") 

pp16_loci <- readLines("../input_data/pp16_loci.txt")

trios_pp16_cpi <- trios_df %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    summarise(n_loci = n(), 
	      n_exclusions = sum(pi == 0),
	      cpi = prod(pi_adj)) %>%
    ungroup() %>%
    mutate(marker_set = "pp16")

trios_sets <- 
    bind_rows(trios_codis_cpi, trios_ident_cpi, trios_pp16_cpi) %>%
    arrange(case_no, trio, desc(n_loci), marker_set)

trios_exclusions <- trios_sets %>%
    filter(n_exclusions > 2)

trios_inconclusive <- trios_sets %>%
    filter(n_exclusions < 3, cpi < 10000)

trios_inclusions <- trios_sets %>%
    filter(n_exclusions < 3, cpi >= 10000)


# total trios
total_trios <- 
    left_join(trios_sets %>% count(marker_set),
	      trios_exclusions %>% count(marker_set), by = c("marker_set")) %>%
    rename(total = n.x, exclusion = n.y)

write_tsv(total_trios, "./total_trios.tsv")


# inclusion to exclusion
left_join(trios_exclusions, trios_cpi, by = c("case_no", "trio")) %>%
    filter(n_exclusions.y < 3, cpi.y >= 10000) %>%
    count(marker_set.x)

# exclusion to inclusion
exc_to_inc <- 
    left_join(trios_inclusions, trios_cpi, by = c("case_no", "trio")) %>%
    filter(n_exclusions.y > 2) %>%
    count(marker_set.x) %>%
    mutate(before = "exclusion", after = "inclusion")

# exclusion to inconclusive
exc_to_inconc <- 
    left_join(trios_inconclusive, trios_cpi, by = c("case_no", "trio")) %>%
    filter(n_exclusions.y > 2) %>%
    count(marker_set.x) %>%
    mutate(before = "exclusion", after = "inconclusive")

# inclusion to inconclusive
inc_to_inconc <- 
    left_join(trios_inconclusive, trios_cpi, by = c("case_no", "trio")) %>%
    filter(n_exclusions.y < 3, cpi.y >= 10000) %>%
    count(marker_set.x) %>%
    mutate(before = "inclusion", after = "inconclusive")

reduction_effect_df <- bind_rows(exc_to_inc, exc_to_inconc, inc_to_inconc) %>%
    select(before, after, marker_set = marker_set.x, n) %>%
    complete(before, after, marker_set, fill = list(n = 0)) %>%
    filter(before != after)

write_tsv(reduction_effect_df, "./reduction_of_loci.tsv")


# separate sets
trios_exclusion_codis <- trios_df %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18, sum(pi == 0) > 2) %>%
    ungroup() %>%
    select(-pi, -pi_adj)

trios_exclusion_ident <- trios_df %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(pi == 0) > 2) %>%
    ungroup() %>%
    select(-pi, -pi_adj)

trios_exclusion_pp16 <- trios_df %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(pi == 0) > 2) %>%
    ungroup() %>%
    select(-pi, -pi_adj)

write_tsv(trios_df, "./trios_pis.tsv") 
write_tsv(trios_exclusion_codis, "./trios_exclusion_codis.tsv") 
write_tsv(trios_exclusion_ident, "./trios_exclusion_ident.tsv") 
write_tsv(trios_exclusion_pp16, "./trios_exclusion_pp16.tsv") 

