sformat <- function(string, format){
  return(do.call(glue::glue, c(string, format)))
}