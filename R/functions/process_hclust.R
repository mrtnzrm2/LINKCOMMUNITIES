process.hclust <- function(net, net.cluster,  nodes, type = 'directed'){
  
  height <- net.cluster$height[!duplicated(net.cluster$height)]
  leaves <- length(height)
  
  Dc.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  Nmin.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  Nmax.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  K.view <- matrix(0, nrow = leaves + 1, ncol = 1)
  NEC.view <- matrix(0, nrow = leaves + 1 , ncol = 1)
  NAC.view <- matrix(nodes, nrow = leaves + 1, ncol = 1)
  
  K.view[1] <- max(net$id)
  
  for (i in 1:leaves){
    
    cluster.h <- cutree(net.cluster, h=height[i])
    clusters <- unique(cluster.h)
    K.view[i+1] <- length(clusters)
    Dc <- matrix(0, nrow = length(clusters), ncol = 1)
    Mc <- matrix(0, nrow = length(clusters), ncol = 1)
    Nc <- matrix(2, nrow = length(clusters), ncol = 1)
    NEC <- matrix(0, nrow = length(clusters), ncol = 1)
    
    for (ii in 1:length(clusters)){
      net.cls <- net[which(cluster.h == clusters[ii]),c('source', 'target')]
      Mc[ii] <- nrow(net.cls)
      sor <- unique(net.cls$source)
      tar <- unique(net.cls$target)
      ndc <- length(unique(c(sor,tar)))
      
      if (ndc > 2 && Mc[ii] > 1){
        if (type == 'directed'){
          Dc[ii] <- (Mc[ii]-ndc+1)/((ndc-1)**2)
          NEC[ii] <- 1
          nds <- unique(intersect(sor, tar))
          if (length(nds) > 1){
            NAC.view[i+1] <- NAC.view[i+1] - length(nds) + 1 
          }
          
          Nc[ii] <- ndc
          
        } else if (type == 'undirected'){
          Dc[ii] <- 2*(Mc[ii]-(ndc-1))/((ndc-1)*(ndc-2))
          NEC[ii] <- 1
          nds <- unique(intersect(sor, tar))
          if (length(nds) > 1){
            NAC.view[i+1] <- NAC.view[i+1] - length(nds) + 1 
          }
          
          Nc[ii] <- ndc
          
        }
      }
      
    }
    Dc.view[i+1] <- t(Mc)%*%Dc/leaves 
    Nmin.view[i+1] <- min(Nc, na.rm = T)
    Nmax.view[i+1] <- max(Nc, na.rm = T)
    NEC.view[i+1] <- sum(NEC == 1)
  }
  
  height <- c(0, height)
  
  D.network <- data.frame('height' = height,
                          'Dc' = Dc.view,
                          'Nmin' = Nmin.view,
                          'Nmax' = Nmax.view,
                          'K' = K.view,
                          'NEC' = NEC.view,
                          'NAC' = NAC.view)
  return(D.network)
}