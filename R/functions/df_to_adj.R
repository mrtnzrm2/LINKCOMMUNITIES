df.to.adj <- function(net){
  
  M <- max(net$source)
  N <- max(net$target)
  
  A <- matrix(0, nrow = M, ncol = N)
  
  for (i in 1:nrow(net))
    A[net$source[i], net$target[i]] <- net$weight[i]
  
  return(A)
}