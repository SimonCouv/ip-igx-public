call_future <- function(){
  x <- future_map_int(1:10, function(x)Sys.getpid())
  message("completed dummy future call")
  message("PIDs", x)
  x
}