library(tidyverse)

calc_pi <- function(marker, a1_F, a2_F, a1_SP, a2_SP, freqs) {

    ch <- c(a1_F, a2_F)
    af <- c(a1_SP, a2_SP)

    af_hom <- length(unique(af)) == 1
    ch_hom <- length(unique(ch)) == 1

    shared_af_ch <- intersect(af, ch)

    if (length(shared_af_ch)==0) {
	0
    
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


freqs <- read_tsv("../../input_data/allele_frequency.tsv")

# CODIS
duos_codis <- read_tsv("../trios/trios_exclusion_codis.tsv") %>% 
    select(-a1_M, -a2_M)

duos_codis$pi <- duos_codis %>% 
    select(-case_no, -trio) %>% 
    pmap_dbl(calc_pi, freqs = freqs)

duos_codis <- mutate(duos_codis, adj_pi = ifelse(pi == 0, 0.001, pi))

duos_codis_inc <- duos_codis %>%
    group_by(case_no, trio) %>%
    summarise(exclusions = sum(pi == 0),
	      cpi = prod(adj_pi)) %>%
    ungroup() %>%
    filter(cpi >= 10000)


# Ident
duos_ident <- read_tsv("../trios/trios_exclusion_ident.tsv") %>% 
    select(-a1_M, -a2_M)

duos_ident$pi <- duos_ident %>% 
    select(-case_no, -trio) %>% 
    pmap_dbl(calc_pi, freqs = freqs)

duos_ident <- mutate(duos_ident, adj_pi = ifelse(pi == 0, 0.001, pi))

duos_ident_inc <- duos_ident %>%
    group_by(case_no, trio) %>%
    summarise(exclusions = sum(pi == 0),
	      cpi = prod(adj_pi)) %>%
    ungroup() %>%
    filter(cpi >= 10000)


# PP16
duos_pp16 <- read_tsv("../trios/trios_exclusion_pp16.tsv") %>% 
    select(-a1_M, -a2_M)

duos_pp16$pi <- duos_pp16 %>% 
    select(-case_no, -trio) %>% 
    pmap_dbl(calc_pi, freqs = freqs)

duos_pp16 <- mutate(duos_pp16, adj_pi = ifelse(pi == 0, 0.001, pi))

duos_pp16_inc <- duos_pp16 %>%
    group_by(case_no, trio) %>%
    summarise(exclusions = sum(pi == 0),
	      cpi = prod(adj_pi)) %>%
    ungroup() %>%
    filter(cpi >= 10000)


inclusion_df <- bind_rows(codis = duos_codis_inc, 
			  identifiler = duos_ident_inc,
			  pp16 = duos_pp16_inc, .id = "marker_set")

write_tsv(inclusion_df, "./false_inclusions.tsv")

