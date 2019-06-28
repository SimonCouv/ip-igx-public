post_process_drake <- function(..., plan) {
  # based on https://ropenscilabs.github.io/drake-manual/plans.html#the-types-of-transformations
  
  args <- list(...)
  names(args) <- all.vars(substitute(list(...)))
  trace <- filter(plan, target %in% names(args))
  # Do post-processing with args and trace.
}