# see demo(error.catching) http://r.789695.n4.nabble.com/capturing-warnings-within-loops-so-I-know-the-iterations-where-warnings-occurred-td4676786.html
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  M <- NULL
  m.handler <- function(m){ # warning handler
    M <<- m
    # invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler,
                                   message = m.handler),
       warning = W,
       message = M)
}