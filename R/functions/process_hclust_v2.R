process.hclust.v2 <- function(path, net, net.cluster,  nodes, labels, type = 'directed'){
  
  source(sprintf('%s/R/functions/parallelization_height.R', path))
  
  y <- parallelization.height(net.cluster$height)
  height <- y$height
  k.step <- y$k
  leaves <- length(height)
  
  Dc.view <- matrix(0, nrow = leaves+1, ncol = 1)
  Nmin.view <- matrix(0, nrow = leaves+1, ncol = 1)
  Nmax.view <- matrix(0, nrow = leaves+1, ncol = 1)
  K.view <- matrix(0, nrow = leaves+1, ncol = 1)
  NEC.view <- matrix(0, nrow = leaves+1, ncol = 1)
  NAC.view <- matrix(nodes, nrow = leaves+1, ncol = 1)
  
  K.view[1] <- max(net$id[net$source <= max(net$target)])
  
  for (i in 1:leaves){
    
    cluster.k <- cutree(net.cluster, k=k.step[i])
    clusters <- unique(cluster.k)
    
    K.view[i+1] <- k.step[i]
    
    Dc <- matrix(0, nrow = k.step[i], ncol = 1)
    Mc <- matrix(0, nrow = k.step[i], ncol = 1)
    Nc <- matrix(2, nrow = k.step[i], ncol = 1)
    NEC <- matrix(0, nrow = k.step[i], ncol = 1)
    
    # if (k.step[i] <= 282)
    #   print(sprintf('k=%i', k.step[i]))
    
    for (ii in 1:k.step[i]){
      net.cls <- net[which(cluster.k == clusters[ii]),c('source', 'target')]
      Mc[ii] <- nrow(net.cls)
      sor <- unique(net.cls$source)
      tar <- unique(net.cls$target)
      ndc <- length(unique(c(sor,tar)))
      
      if (ndc > 2 && Mc[ii] > 1){
        
        if (type == 'directed'){
          Dc[ii] <- (Mc[ii]-ndc+1)/((ndc-1)**2)
          
          nds <- unique(intersect(sor, tar))
          if (length(nds) > 1){
            NEC[ii] <- 1
          }
          
          nd.insct <- intersect(sor,tar)
          nd <- nd.insct %>% length()
          if (length(nd.insct) > 1 && Mc[ii] >= 0.8*nd*(nd-1)){
            NAC.view[i+1] <- NAC.view[i+1]  + 1 - length(nd.insct)
            # if (k.step[i] <= 282){
            #   print(clusters[ii])
            #   print(labels[nd.insct])
            # }
          }
          
          Nc[ii] <- ndc
          
        } else if (type == 'undirected'){
          Dc[ii] <- 2*(Mc[ii]-(ndc-1))/((ndc-1)*(ndc-2))
          
          nds <- unique(intersect(sor, tar))
          if (length(nds) > 1){
            NEC[ii] <- 1
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