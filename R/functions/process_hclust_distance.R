calculate.s.s <- function(s, M, n=2){
  s.s <- 0
  ns <- length(s)
  if (ns > 1){
    if (ns > n){
      for (i in 2:n){
        s.s <- s.s + (s[i]/s[i-1])*((s[i]+s[i-1])/M)
      }
      s.s <- s.s/n
    } else{
      for (i in 2:ns){
        s.s <- s.s + (s[i]/s[i-1])*((s[i]+s[i-1])/M)
      }
      s.s <- s.s/ns
    }
    
  } else
    s.s <- NA
  return(s.s)
}

process.hclust.dist <- function(net.cluster){
  M <- net.cluster$height %>% length()
  M <- M + 1
  height <- net.cluster$height[!duplicated(net.cluster$height)]
  leaves <- length(height)
  s.s.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  K.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  
  for (i in 1:leaves){
    cluster.h <- cutree(net.cluster, h=height[i])
    clusters <- unique(cluster.h)
    K.view[i+1] <- length(clusters)
    s.s <- cluster.h %>% table() %>% sort(decreasing = T)
    s.s <- s.s %>% unname()
    s.s.view[i+1] <- calculate.s.s(s.s, M, n=3)
  }
  
  height <- c(0, height)
  
  D.network <- data.frame('height' = height,
                          "SS"= s.s.view,
                          'K' = K.view)
  return(D.network)
}