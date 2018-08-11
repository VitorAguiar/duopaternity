library(tidyverse)

calc_pi <- function(marker, a1_F, a2_F, a1_SP, a2_SP, freqs) {

    ch <- c(a1_F, a2_F)
    af <- c(a1_SP, a2_SP)

    af_hom <- length(unique(af)) == 1
    ch_hom <- length(unique(ch)) == 1

    shared_af_ch <- intersect(af, ch)

    if (length(shared_af_ch)==0) {
	0.001
    
    } else if (af_hom && ch_hom && all(af == ch)) {
	a <- unique(af)
	fa <- freqs$f[freqs$marker == marker & freqs$allele == a]
	1/fa
    
    } else if (((!af_hom && ch_hom) || (af_hom && !ch_hom)) && length(shared_af_ch)==1) {
	a <- shared_af_ch
	fa <- freqs$f[freqs$marker == marker & freqs$allele == a]
	1/(2*fa)

    } else if (!af_hom && !ch_hom && length(shared_af_ch)==2) {
	a <- af[1]
	b <- af[2]
	fa <- freqs$f[freqs$marker == marker & freqs$allele == a]
	fb <- freqs$f[freqs$marker == marker & freqs$allele == b]
	(fa + fb)/(4*fa*fb)

    } else if (!af_hom && !ch_hom && length(shared_af_ch)==1) {
	a <- shared_af_ch
	fa <- freqs$f[freqs$marker == marker & freqs$allele == a]
	1/(4*fa)

    } else {
	stop("case not predicted")
    }
}


# input data

freqs <- read_tsv("../input_data/allele_frequency.tsv")

# All trios
duos <- read_tsv("../trios/exclusion_trios.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)

duos$pi <- duos %>% select(-case_no, -trio) %>% pmap_dbl(calc_pi, freqs = freqs)

# Filter for Codis 
codis_loci <- readLines("../input_data/codis_loci.txt")

inclusions_codis <- duos %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

count(duos, case_no, trio)
count(inclusions, case_no, trio)


# Filter for Identifiler 
identifiler_loci <- readLines("../input_data/identifiler_loci.txt")

inclusions_identifiler <- duos %>%
    filter(marker %in% identifiler_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

# Filter for pp16 
pp16_loci <- readLines("../input_data/pp16_loci.txt")

inclusions_pp16 <- duos %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup()


# Preselected exclusion trios based on CODIS
duos_codis <- read_tsv("../trios/exclusion_trios_codis.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)

duos_codis$pi <- 
    duos_codis %>% select(-case_no, -trio) %>% pmap_dbl(calc_pi, freqs = freqs)

inclusions_codis_v2 <- duos_codis %>%
    group_by(case_no, trio) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup()


# Preselected exclusion trios based on Identifiler
duos_identifiler <- read_tsv("../trios/exclusion_trios_identifiler.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)

duos_identifiler$pi <- 
    duos_identifiler %>% select(-case_no, -trio) %>% pmap_dbl(calc_pi, freqs = freqs)

inclusions_identifiler_v2 <- duos_identifiler %>%
    group_by(case_no, trio) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup()

# Preselected exclusion trios based on PP16
duos_pp16 <- read_tsv("../trios/exclusion_trios_pp16.tsv") %>% 
    select(-a1_M, -a2_M, -pi, -cpi)

duos_pp16$pi <- 
    duos_pp16 %>% select(-case_no, -trio) %>% pmap_dbl(calc_pi, freqs = freqs)

inclusions_pp16_v2 <- duos_pp16 %>%
    group_by(case_no, trio) %>%
    mutate(cpi = prod(pi)) %>%
    filter(cpi >= 10000) %>%
    ungroup():wq


