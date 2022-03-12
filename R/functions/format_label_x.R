format.label.x <- function(labels){
  for (i in 1:length(labels)){
    if (grepl("X", labels[i], fixed = T))
      labels[i] <- stringr::str_replace(labels[i], "X", "")
  }
  return(labels)
}