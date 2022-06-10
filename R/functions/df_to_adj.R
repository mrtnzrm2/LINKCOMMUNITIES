df.to.adj <- function(net) {
  M <- max(net$source)
  N <- max(net$target)
  A <- matrix(0, nrow = M, ncol = N)
  for (i in seq_len(nrow(net)))
    A[net$source[i], net$target[i]] <- net$weight[i]
  return(A)
}