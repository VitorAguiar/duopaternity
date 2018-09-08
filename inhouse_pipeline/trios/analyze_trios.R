library(tidyverse)

calc_pi <- function(marker, a1_M, a2_M, a1_F, a2_F, a1_SP, a2_SP, freqs) {

  mo <- c(a1_M, a2_M)
  ch <- c(a1_F, a2_F)
  af <- c(a1_SP, a2_SP)
  
  af_hom <- length(unique(af)) == 1L
  ch_hom <- length(unique(ch)) == 1L
  mo_hom <- length(unique(mo)) == 1L

  shared_af_ch <- intersect(af, ch)
  shared_af_mo <- intersect(af, mo)
  shared_mo_ch <- intersect(mo, ch)

  if (!ch_hom && length(shared_mo_ch)==1 && all(shared_mo_ch == shared_af_ch)) {
      return(0)
  }

  if (length(shared_af_ch) == 0) {
      return(0)
  }

  if (af_hom && ch_hom && all(af==ch) && all(af %in% mo)) {
      a <- unique(af)
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      1/fa

  } else if (af_hom && !ch_hom && length(shared_af_ch)==1 && length(shared_af_mo)==0) {
      a <- unique(af)
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      1/fa

  } else if (!af_hom && ch_hom && all(shared_af_ch %in% mo)) {
      a <- unique(ch)
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      1/(2*fa)

  } else if (!af_hom && !ch_hom && length(shared_af_ch)==2 && length(shared_mo_ch)==1) {
      a <- ch[!ch %in% mo & ch %in% af]
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      1/(2*fa)

  } else if (!af_hom && !ch_hom && length(shared_af_ch)==1 && all(!shared_af_ch %in% mo))  {
      a <- ch[!ch %in% mo & ch %in% af]
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      1/(2*fa)

  } else if (!mo_hom && ! ch_hom && length(shared_mo_ch)==2 && all(af %in% mo)) {
      a <- ch[1]
      b <- ch[2]
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      fb <- freqs$f[freqs$marker == marker & freqs$allele == b] 
      1/(fa + fb)

  } else if (!mo_hom && !ch_hom && !af_hom && length(shared_mo_ch)==2 && length(shared_af_ch)==1) {
      a <- ch[1]
      b <- ch[2]
      fa <- freqs$f[freqs$marker == marker & freqs$allele == a] 
      fb <- freqs$f[freqs$marker == marker & freqs$allele == b] 
      1/(2*(fa + fb))
  
  } else {
    stop("case not predicted")
  }
}
 
trios <- read_tsv("../../input_data/integrated_data.tsv")
names(trios) <- sub("allele\\.", "a", names(trios))

freqs <- read_tsv("../../input_data/allele_frequency.tsv")

trios$pi <- trios %>% 
    select(-case_no, -trio) %>% 
    pmap_dbl(calc_pi, freqs = freqs)

trios <- mutate(trios, adj_pi = ifelse(pi == 0, 0.001, pi))


codis_loci <- sort(readLines("../../input_data/codis_loci.txt"))

trios_exclusion_codis <- trios %>%
    filter(marker %in% codis_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 18, sum(pi == 0) > 2) %>% 
    ungroup() %>%
    select(-pi, -adj_pi)

ident_loci <- sort(readLines("../../input_data/identifiler_loci.txt"))

trios_exclusion_ident <- trios %>%
    filter(marker %in% ident_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(pi == 0) > 2) %>% 
    ungroup() %>%
    select(-pi, -adj_pi)

pp16_loci <- sort(readLines("../../input_data/pp16_loci.txt"))

trios_exclusion_pp16 <- trios %>%
    filter(marker %in% pp16_loci) %>%
    group_by(case_no, trio) %>%
    filter(n() == 15, sum(pi == 0) > 2) %>% 
    ungroup() %>%
    select(-pi, -adj_pi)

write_tsv(trios, "./trios_pi.tsv")
write_tsv(trios_exclusion_codis, "./trios_exclusion_codis.tsv")
write_tsv(trios_exclusion_ident, "./trios_exclusion_ident.tsv")
write_tsv(trios_exclusion_pp16, "./trios_exclusion_pp16.tsv")

