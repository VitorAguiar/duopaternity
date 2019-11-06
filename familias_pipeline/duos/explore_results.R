library(tidyverse)

trios <- read_tsv("../../input_data/integrated_data.tsv")

shared16 <- trios %>% filter(marker == "D12S391", (ch_1 == 16.3 | ch_2 == 16.3) & (af_1 == 16.3 | af_2 == 16.3))
notd7 <- trios %>% filter(marker == "D7S820", ch_1 != af_1 & ch_2 != af_1 & ch_2 != af_1 & ch_2 != af_2)

shared16 %>% filter(case_no %in% notd7$case_no)
read_tsv("../trios/trios_cpi.tsv") %>% filter(case_no == 250992)

trios_cpi <- read_tsv("../trios/trios_cpi.tsv") %>% 
    filter(n_loci >= 18, n_exclusions >= 4)

duos_pi_r001 <- read_tsv("./duos_pi_r001.tsv")
duos_cpi_r001 <- read_tsv("./duos_cpi_r001.tsv")
duos_inclusion_r001 <- read_tsv("./duos_inclusion_r001.tsv", col_types = "dcdd")

duos_pi_sw <- read_tsv("./duos_pi_stewise.tsv")
duos_cpi_sw <- read_tsv("./duos_pi_stewise.tsv")
duos_inclusion_sw <- read_tsv("./duos_inclusion_stepwise.tsv")

inclusion_df <- list(r001 = duos_inclusion_r001, stepwise = duos_inclusion_sw) %>%
    bind_rows(.id = "model") %>%
    left_join(trios_cpi, by = c("case_no", "trio"), suffix = c(".duo", ".trio"))

duos_cpi_r001 %>% filter(case_no %in% inclusion_df$case_no)

x <- duos_pi_r001 %>% filter(case_no %in% inclusion_df$case_no)
y <- duos_pi_sw %>% filter(case_no %in% inclusion_df$case_no)

left_join(x, y, by = c("case_no", "trio", "marker"), suffix = c(".r001", ".sw")) %>%
    print(n = Inf)
