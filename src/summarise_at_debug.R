library(dplyr)
library(drake)

test <- tibble(abc=c(NA, 4,5 ), foo=1:3)

test %>% dplyr::summarise_at(.vars=c("abc", "foo"), .funs=mean)   # always works

test %>% dplyr::summarise_at(.vars=vars(-one_of("abc")), .funs=mean)   # works after restarting

test %>% dplyr::summarise_at(.vars=vars(starts_with("foo")), .funs=mean)  # broken in original interactive session, works in reprex and after restarting session

test %>% dplyr::select(starts_with("foo"))   # always works

session_info <- sessionInfo()
saveRDS(session_info, "summarise_at_error_info.RDS")
