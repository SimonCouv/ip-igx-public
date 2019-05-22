# parse file name
args <- commandArgs(trailingOnly = TRUE)
fp <- args[1]
stopifnot(grepl(tolower(fp), pattern = "\\.rmd$"))
fn <- fs::path_ext_remove(fp)
script <- paste0(fn, ".R")

# purl
knitr::purl(fp)

# source
source(script)

# deleted script
fs::file_delete(script)
