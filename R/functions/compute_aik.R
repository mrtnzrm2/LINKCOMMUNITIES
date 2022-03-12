compute.aik <- function(net, N, mode="ALPHA"){
  source("functions/aki.R")
  source("functions/aik.R")
  source("functions/jaccard_p_fast.R")
  
  aik.matrix <- aik(net, mode=mode)
  AIK <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N){
    for (j in 1:N){
      if (i < j){
        AIK[i,j] <- jaccard.p.fast(aik.matrix[i,], aik.matrix[j,])
        
      }
    }
  }
  AIK <- AIK + t(AIK)
  return(AIK)
}