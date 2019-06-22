library(tidyverse)


falseinc_r001 <- read_tsv("./false_inclusions_r001.tsv")
falseinc_strb <- read_tsv("./false_inclusions_strbase.tsv")

falseinc_df <- list(r001 = falseinc_r001, strbase = falseinc_strb) %>%
    bind_rows(.id = "model")

trios <- read_tsv("../trios/trios_pis.tsv")

# ident
duos_ident_r001 <- read_tsv("./duos_ident_pi_r001.tsv")
duos_ident_strb <- read_tsv("./duos_ident_pi_strbase.tsv")
ident_loci <- sort(readLines("../../input_data/identifiler_loci.txt"))

trios %>%
    filter(case_no == 253626, trio == "M1_F1_SP1")

trios %>%
    filter(case_no == 253626, trio == "M1_F1_SP1", marker %in% ident_loci)

duos_ident_r001 %>%
    filter(case_no == 253626, trio == "M1_F1_SP1")

duos_ident_strb %>%
    filter(case_no == 253626, trio == "M1_F1_SP1")

# pp16
duos_pp16_r001 <- read_tsv("./duos_pp16_pi_r001.tsv")
duos_pp16_strb <- read_tsv("./duos_pp16_pi_strbase.tsv")
pp16_loci <- sort(readLines("../../input_data/pp16_loci.txt"))

trios %>%
    filter(case_no == 117770, trio == "M1_F1_SP1", marker %in% pp16_loci)

duos_pp16_r001 %>%
    filter(case_no == 117770, trio == "M1_F1_SP1")

duos_pp16_strb %>%
    filter(case_no == 117770, trio == "M1_F1_SP1")

#
duos_pp16_r001 %>%
    filter(case_no == 155180, trio == "M1_F1_SP1")

duos_pp16_strb %>% 
    filter(case_no == 155180, trio == "M1_F1_SP1") %>%
    summarise(cpi = prod(pi),
	      n_exclusions = sum(exclusion))

#
duos_codis_r001 <- read_tsv("./duos_codis_pi_r001.tsv")
codis_loci <- sort(readLines("../../input_data/codis_loci.txt"))

trios %>%
    filter(case_no == 250992, trio == "M1_F1_SP1", marker %in% codis_loci) %>%
    summarise(n_exclusions = sum(exclusion))

duos_codis_r001 %>%
    filter(case_no == 250992, trio == "M1_F1_SP1") %>%
    summarise(n_exclusions = sum(exclusion))

#
duos_codisplus_r001 <- read_tsv("./duos_codisplus_pi_r001.tsv")
codisplus_loci <- sort(readLines("../../input_data/codisplus_loci.txt"))

trios %>%
    filter(case_no == 250992, trio == "M1_F1_SP1", marker %in% codisplus_loci) %>%
    summarise(n_exclusions = sum(exclusion))

duos_codisplus_r001 %>%
    filter(case_no == 250992, trio == "M1_F1_SP1") %>%
    summarise(n_exclusions = sum(exclusion))

trios %>%
    filter(case_no == 250992, trio == "M1_F1_SP1")


