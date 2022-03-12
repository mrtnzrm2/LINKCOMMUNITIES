format.labels <- function(labels){
  for (i in 1:length(labels)){
    if (grepl("/",labels[i], fixed = T)){
      labels[i] <- stringr::str_replace(labels[i], "/", ".")
    }
  }
  return(labels)
}