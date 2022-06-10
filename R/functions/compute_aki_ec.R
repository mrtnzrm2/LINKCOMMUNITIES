compute.aki.ec <- function(net, N, mode="ALPHA") {
  source("functions/aki_ec.R")
  source("functions/jaccard_p_fast.R")
  aki.matrix <- aki_ec(net, mode = mode)
  AKI <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j) {
        AKI[i, j] <- jaccard.p.fast(aki.matrix[i, ], aki.matrix[j, ])
      }
    }
  }
  AKI <- AKI + t(AKI)
  return(AKI)
}