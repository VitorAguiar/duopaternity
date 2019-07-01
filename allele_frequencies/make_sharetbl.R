suppressPackageStartupMessages(library(tidyverse))

opts <- commandArgs(TRUE)
profiles_df <- opts[1]
chunk       <- as.integer(opts[2])
nchunks     <- as.integer(opts[3])
outPrefix   <- opts[4]

strtbl <- suppressMessages(read_tsv(profiles_df, progress = FALSE)) %>%
    mutate(ix = group_indices(., id)) %>%
    select(ix, everything()) %>%
    arrange(ix)

indices <- tibble(index = unique(strtbl$ix)) %>%
    mutate(ck = ntile(index, nchunks))

indices_chunk <- filter(indices, ck == chunk)

out_list <- vector("list", nrow(indices_chunk))
counter <- 1L

for (i in indices_chunk$index) {

    dat_i <- filter(strtbl, ix == i) %>% 
	select(-ix)
    
    dat_iplus <- filter(strtbl, ix > i) %>% 
	select(-ix)
    
    str_i <- left_join(dat_i, dat_iplus, by = "marker") %>% 
	filter(id.x != id.y) %>% 
	drop_na()

    out_list[[counter]] <- str_i %>%
	mutate(i1 = as.integer(hap.x == hap.y & allele.x == allele.y),
	       i2 = as.integer(hap.x != hap.y & allele.x == allele.y)) %>%
	group_by(id.x, id.y, marker) %>%
	summarise(i1 = sum(i1), i2 = sum(i2)) %>%
	ungroup() %>%
	mutate(shr = pmax(i1, i2)) %>%
	select(-i1, -i2) %>%
	group_by(id.x, id.y) %>%
	summarise(total = sum(shr == 2),
		  partial = sum(shr == 1)) %>%
	ungroup()

    counter <- counter + 1L
}

out_tbl <- bind_rows(out_list) %>%
    filter(total > 13)

out <- sprintf("%s_%d.tsv", outPrefix, chunk)
write_tsv(out_tbl, out)
