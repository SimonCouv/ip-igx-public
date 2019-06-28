library(stringr)
library(BiocManager)

ncpu <- commandArgs(trailingOnly = TRUE)[1]
target <- commandArgs(trailingOnly = TRUE)[2]
# target <- "src/imports_slurm.R"
if (is.na(ncpu))
  ncpu <- 1

# function definition ------------------------------------------
install_pkgs_from_library_calls <- function(library_calls_fp,
                                            lib_path = .libPaths()[1],
                                            ncpu=1){
  # determine pkgs to install
  lines <- readLines(library_calls_fp)
  library_pkgs <- str_match(lines, pattern = "^library\\(([^,]+).*\\)$")[,2]
  to_install <-setdiff(library_pkgs, installed.packages()[,'Package'])
  
  # install, if required
  if (length(to_install)>0){
    message("installing ", length(to_install), " packages")
    install.packages(to_install, 
                     repos = BiocManager::repositories(),   # includes CRAN,
                     Ncpus = ncpu
    )
    message("installation complete")
  } else
    message("all required packages are already installed. Exiting.")
  
}

# Execution ------------------------------------------
# call function
install_pkgs_from_library_calls(library_calls_fp = target,
                                ncpu=ncpu)