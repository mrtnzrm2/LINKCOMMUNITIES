assign.commship <- function(net, net.cluster, k,  bh, kactive=T){
  if (kactive){
    net$commship <- cutree(net.cluster, k=k)
  } else{
    net$commship <- cutree(net.cluster, h=bh)
  }
  return(net)
}