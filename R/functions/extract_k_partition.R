extract.k.partition <- function(net, net.cluster, nodes, k){

  source('functions/linkcomm_parameters.R')
  
  net$commship <-  cutree(net.cluster, k=k)
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  net <- net[which(net$commship != -1),]
  
  node.historia <- rep(-1, nodes)
  
  coms <- sort(unique(net$commship))
  
  for (com in coms){
    net.com <- net[which(net$commship == com),]
    nds <- sort(unique(intersect(net.com$source, net.com$target)))
    if (length(nds) > 1)
      node.historia[nds] <- com
  }
  return(node.historia)
}