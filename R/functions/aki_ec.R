aki_ec <- function(net, self.loop = T, mode="") {
  net <- net[net$source <= 50, ]
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  matrix <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  for (i in 1:nodes.t) {
    matrix[i, net$source[which(net$target == i)]] <-
      net$weight[which(net$target == i)]
    if (self.loop) {
      if (pracma::strcmp(mode, "BETA")) {
        matrix[i, i] <-
          mean(net.source$weight[which(net.source$source == i)], na.rm = T)
      } else if (pracma::strcmp(mode, "ALPHA")) {
        matrix[i, i] <- mean(net$weight[which(net$target == i)], na.rm = T)
      }
    }
  }
  return(matrix)
}