library(tidyverse)

argv <- commandArgs(TRUE)
groups_file <- argv[1]
chunk <- as.integer(argv[2])
nchunks <- as.integer(argv[3])
outPrefix <- argv[4]

profile_groups <- read_tsv(groups_file)

n_inds <- n_distinct(profile_groups$id)

indices <- tibble(id = sort(unique(profile_groups$id)),
		  ck = ntile(id, nchunks))

indices_chunk <- filter(indices, ck == chunk)

groups_chunk <- profile_groups %>%
    group_by(gid, marker) %>%
    filter(any(id %in% indices_chunk$id)) %>%
    ungroup() %>%
    arrange(gid, id, marker)

share_list <- vector("list", nrow(indices_chunk))

for (i in indices_chunk$id) {

    ind_i <- groups_chunk %>%
	group_by(gid, marker) %>%
	filter(any(id == i)) %>%
	ungroup() %>%
	mutate(marker = sub("\\.[12]$", "", marker))

    share_list[[i]] <- ind_i %>% 
	filter(id != i) %>% 
	count(id, marker) %>%
	group_by(id) %>%
	summarise(total = sum(n == 2L),
		  partial = sum(n == 1L)) %>%
	ungroup() %>%
	filter(total >= 13)
}

out_tbl <- bind_rows(share_list, .id = "id1") %>%
    rename(id2 = id)

outname <- sprintf("%s_%d.tsv", outPrefix, chunk)

write_tsv(out_tbl, outname)
