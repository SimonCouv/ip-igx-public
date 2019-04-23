clean_rmd_for_render <- function(rmd, out_name, comment_regex){
  rmd_f <- file(rmd, open = "r")
  out_f <- file(out_name, open = "w")
  while (TRUE) {
    line <- readLines(rmd_f, n=1)
    
    if(length(line)==0)
      break
    
    if (any(str_detect(line, pattern = comment_regex))){
      line <- paste0("# ", line)
    }
    writeLines(line, con = out_f)

  }
  return(out_name)
}