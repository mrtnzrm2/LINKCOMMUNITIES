find.height <- function(Dnetwork, curve='Dview', max.wire = 1, min.wire = 0, 
                        l_continuous=T, jump=0.07){
  
  f <- Dnetwork[,curve]
  best.f <- -Inf
  best.height <- 0
  
  max.h <- max(f)
  min.h <- min(f)
  
  if (max.h < abs(min.h))
    max.h <- abs(min.h)
  
  if (l_continuous){
    for (i in 2:length(Dnetwork$height)){
      if (i < length(Dnetwork$height) && (f[i]-f[i-1])/max.h < jump){
        if ((f[i] >= best.f && i/length(Dnetwork$height) <= max.wire && i/length(Dnetwork$height) >= min.wire)){
          best.height <- Dnetwork$height[i]
          best.f <- f[i]
        }
      } else if (i < length(Dnetwork$height) && (f[i]-f[i-1])/max.h >= jump){
        break
      }
    }
  }else{
    for (i in 2:length(Dnetwork$height)){
      if (i < length(Dnetwork$height)){
        if ((f[i] >= best.f && i/length(Dnetwork$height) <= max.wire && i/length(Dnetwork$height) >= min.wire)){
          best.height <- Dnetwork$height[i]
          best.f <- f[i]
        }
      } 
    }
  }
  
  return(best.height)
}