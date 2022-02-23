with.potholes <- function(a){
  
  a.notrees <- a[a != -1]
  
  element.a <- c()
  size.a <- c()
  for(e in 1:length(a)){
    if (!(a[e] %in% element.a) || a[e] == -1){
      element.a <- c(element.a, a[e])
      if (a[e] != -1){
        size.a <- c(size.a, sum(a.notrees == a[e]))
      } else {
        size.a <- c(size.a, 1)
      }
    }
  }
  
  return(list(element.a, size.a))
}