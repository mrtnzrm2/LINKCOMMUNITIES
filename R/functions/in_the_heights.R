in.the.heights <- function(merde, nodes){
  
  height <- matrix(0, nrow = nodes, ncol = 1)
  inv.height <- sort(merde$height, decreasing = T)
  for (e in (nodes-1):1){
    tree.nodes <- cutree(merde, k=e)
    comms <- unname(tree.nodes[duplicated(tree.nodes)])
    tree.nodes <- which(tree.nodes %in% comms)
    height[tree.nodes] <- height[tree.nodes] + inv.height[e] 
    # for (nd in tree.nodes){
    #   if (height[nd] == 0){
    #     height[nd] <- inv.height[e] 
    #   } else{
    #     height[nd] <- height[nd] + (inv.height[e] - height[nd])/(nodes - e + 1)
    #   }
    # }
  }
  height <- height - min(height)
  return(height)
}