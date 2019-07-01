library(tidyverse)



dups <- read_tsv("./sharetable_all.tsv")

ids <- c("133155_af", "158521_af", "134684_m") 

dups %>% 
    filter(id.x %in% ids | id.y %in% ids) %>%
    mutate(i = 1) %>%
    select(id.x, id.y, i) %>%
    spread(id.y, i)


