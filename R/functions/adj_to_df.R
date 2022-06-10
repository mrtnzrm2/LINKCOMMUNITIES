adj.to.df <- function(A) {
  M <- dim(A)[1]
  N <- dim(A)[2]
  weights <- pracma::Reshape(t(A), M * N, 1)
  target <- rep(1:N, N = M)
  source <- rep(1:M, each = N)
  df <- data.frame(source = source, target = target, weight = weights)
  return(df)
}