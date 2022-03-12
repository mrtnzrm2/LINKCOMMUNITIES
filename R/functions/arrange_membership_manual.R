arrange.membership.manual <- function(id){
  for (i in 1:length(id)){
    if (id[i] == 4)
      id[i] <- 2
    else if (id[i] == 3)
      id[i] <- 4
    else if (id[i] == 2)
      id[i] <- 6
    else if (id[i] == 6)
      id[i] <- 3
  }
  return(id)
}