arrange.membership.manual <- function(id, arragement) {
  for (i in seq_len(length(id))) {
    for (name in names(arragement)) {
      if (id[i] == as.numeric(name))
        id[i] <- arragement[name]
      else if (id[i] == arragement[name])
        id[i] <- as.numeric(name)
    }
  }
  return(id)
}