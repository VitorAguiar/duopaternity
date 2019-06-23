library(tidyverse)

#opts <- commandArgs(TRUE)
#profiles_df <- opts[1]
#outPrefix   <- opts[2]

profiles_df <- "./str_parents.tsv"
outPrefix <- "./sharetable"

strtbl <- read_tsv(profiles_df) %>%
    mutate(ix = group_indices(., id)) %>%
    select(ix, everything()) %>%
    arrange(ix)

for (i in unique(strtbl$ix)) {

    dat_i <- filter(strtbl, ix == i) %>%
	left_join(filter(strtbl, ix > i), by = "marker") %>%
	drop_na()

    share_tbl_i <- dat_i %>%
	mutate(i1 = as.integer(hap.x == hap.y & allele.x == allele.y),
	       i2 = as.integer(hap.x != hap.y & allele.x == allele.y)) %>%
	group_by(id.x, id.y, marker) %>%
	summarise(i1 = sum(i1), i2 = sum(i2)) %>%
	ungroup() %>%
	mutate(shr = pmax(i1, i2)) %>%
	select(-i1, -i2) %>%
	group_by(id.x, id.y) %>%
	summarise(total = mean(shr == 2),
		  partial = mean(shr == 1)) %>%
	ungroup()

    out_i <- sprintf("%s_%d.tsv", outPrefix, i)
    write_tsv(share_tbl_i, out_i)
    print(i)
}

