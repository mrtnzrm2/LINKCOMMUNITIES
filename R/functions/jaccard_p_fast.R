jaccard.p.fast <- function(x, y){
  
  n <- length(x)
  
  x.m <- matrix(rep(x, n), ncol = n, byrow = T)
  x.m <- x.m/diag(x.m)
  x.m[is.nan(x.m)] <- Inf
  
  y.m <- matrix(rep(y, n), ncol = n, byrow = T)
  y.m <- y.m/diag(y.m)
  y.m[is.nan(y.m)] <- Inf
  
  return(sum(1/apply(pmax(x.m, y.m), c(1), sum)))
}