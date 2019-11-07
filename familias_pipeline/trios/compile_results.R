library(tidyverse)

trios <- map_df(sprintf("trios_results_batch%d.tsv", 1:8), read_tsv)

write_tsv(trios, "results.tsv")
