with.potatoes <- function(a){
  element.a <- c()
  size.a <- c()
  for(e in 1:length(a)){
    if (!(a[e] %in% element.a)){
      element.a <- c(element.a, a[e])
      size.a <- c(size.a, sum(a == a[e]))
    }
  }
  
  return(list(element.a, size.a))
}