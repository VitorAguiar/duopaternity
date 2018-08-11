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
      return(0.001)
  }

  if (length(shared_af_ch) == 0) {
      return(0.001)
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
 
trios <- read_tsv("../input_data/integrated_data.tsv")
names(trios) <- sub("allele\\.", "a", names(trios))

freqs <- read_tsv("../input_data/allele_frequency.tsv")

trios$pi <- trios %>% 
    select(-case_no, -trio) %>% 
    pmap_dbl(calc_pi, freqs = freqs)


