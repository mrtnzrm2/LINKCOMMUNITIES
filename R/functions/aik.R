aik <- function(net, self.loop = T, mode="") {
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  matrix <- matrix(0, nrow = nodes.s, ncol = nodes.t)
  for (i in 1:nodes.s) {
    matrix[i, net$target[net$source == i]] <- net$weight[net$source == i]
    if (self.loop) {
      if (i <= nodes.t) {
        if (pracma::strcmp(mode, "BETA")) {
          matrix[i, i] <- mean(net$weight[net$target == i], na.rm = T, )
        } else if (pracma::strcmp(mode, "ALPHA")) {
          matrix[i, i] <- mean(net$weight[net$source == i], na.rm = T)
        }
      }
    }
  }
  return(matrix)
}