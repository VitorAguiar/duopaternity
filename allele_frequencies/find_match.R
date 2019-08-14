library(tidyverse)

argv <- commandArgs(TRUE)
groups_file <- argv[1]
chunk <- as.integer(argv[2])
chunk_file <- argv[3]
outPrefix <- argv[4]

find_match <- function(ind, groups_df) 
    groups_df %>%
	group_by(marker, g) %>%
	filter(any(id == ind)) %>%
	ungroup() %>%
	filter(id != ind) %>%
	count(id) %>%
	filter(n > 13)

groups_df <- read_tsv(groups_file)

chunk_df <- read_tsv(chunk_file) %>%
    filter(ck == chunk)

idlist <- as.list(chunk_df$id)
names(idlist) <- unlist(idlist)

out <- map_df(idlist, find_match, groups_df, .id = "id1") %>%
    select(id1, id2 = id, n)

outname <- paste0(outPrefix, "_", chunk, ".tsv")
write_tsv(out, outname)
