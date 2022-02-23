parallelization.height <- function(height){
  steps <- length(height)
  k.step <- steps:1
  pf <- rep(F, steps)
  
  tmp <- height[1]
  for (i in 2:steps){
    if (tmp != height[i]){
      pf[i-1] <- T
      tmp <- height[i]
    }
  }
  return(list(height=height[pf], k=k.step[pf]))
}