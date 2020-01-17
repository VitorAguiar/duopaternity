library(tidyverse)

argv <- commandArgs(TRUE)
groups_file <- argv[1]
chunk <- as.integer(argv[2])
outPrefix <- argv[3]

find_match <- function(ind, str_df) {
    
    nloci <- str_df %>%
        filter(id == ind) %>%
        nrow()
    
    str_df %>%
        filter(any(id == ind)) %>%
        ungroup() %>%
        count(id) %>%
        filter(id != ind) %>%
        filter(n >= (nloci - 1L))
}

groups_df <- read_tsv(groups_file) 

chunk_ids <- filter(groups_df, ck == chunk) %>%
    pull(id) %>%
    unique() %>%
    setNames(., .)

str_df_i <- select(groups_df, id, g) %>%
    group_by(g)

out <- map_df(chunk_ids, find_match, str_df = str_df_i, .id = "ind1")

outname <- paste0(outPrefix, "_", chunk, ".tsv")
write_tsv(out, outname)
