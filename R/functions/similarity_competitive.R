similarity.competitive <- function(net, leaves, self.loop = F, mode=''){
  
  source('functions/jaccard_p_fast.R')
  
  nodes.s <- max(net$source)
  nodes.t <- max(net$target)
  
  aik <- matrix(0, nrow = nodes.t, ncol = nodes.t)
  net.source <- net[which(net$source <= nodes.t),]
  
  for (i in 1:nodes.t){
    aik[i, net.source$target[which(net.source$source == i)]] <- net.source$weight[which(net.source$source == i)]
  
    if (self.loop){
      if (pracma::strcmp(mode, 'BETA')){
        aik[i,i] <- mean(net$weight[which(net$target == i)], na.rm = T,)  
      } else if (pracma::strcmp(mode, 'ALPHA')){
        aik[i,i] <- mean(net.source$weight[which(net.source$source == i)], na.rm = T)
      } 
    }
  }
  
  aki <- matrix(0, nrow = nodes.t, ncol = nodes.s)
  
  for (i in 1:nodes.t){
    aki[i, net$source[which(net$target == i)]] <- net$weight[which(net$target == i)]
    
    if (self.loop){
      if (pracma::strcmp(mode, 'BETA')){
        aki[i,i] <- mean(net.source$weight[which(net.source$source == i)], na.rm = T)
      } else if (pracma::strcmp(mode, 'ALPHA')){
        aki[i,i] <- mean(net$weight[which(net$target == i)], na.rm = T)
      } 
    }
  }
  
  if (pracma::strcmp(mode, 'GAMMA')){
    for (i in 1:nodes.t){
      corr.i <- jaccard.p(aki[i,], aik[i,])
      aki[i,i] <- corr.i
      aik[i,i] <- corr.i
    }
  }
  
  net.sim <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:leaves){
    for (j in 1:leaves){
      if (i < j && net$id[i] <= leaves && net$id[j] <= leaves){
        if (net$source[i] == net$source[j] && net$target[i] != net$target[j]) {
          net.sim[net$id[i],net$id[j]] <- jaccard.p.fast(aki[net$target[i],], aki[net$target[j],])
        } 
        else if (net$source[i] != net$source[j] && net$target[i] == net$target[j]) {
          net.sim[net$id[i],net$id[j]] <- jaccard.p.fast(aik[net$source[i],], aik[net$source[j],])
        }
      }
    }
  }
  net.sim <- net.sim + t(net.sim)
  net.sim[which(net.sim == 0)] <- NA
  return(net.sim)
}