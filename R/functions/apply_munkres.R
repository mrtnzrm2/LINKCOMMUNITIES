apply.munkres <- function(x, y, kmax){
  
  COST <- matrix(0, nrow = kmax, ncol = kmax)
  
  for (i in 1:kmax){
    for (j in 1:kmax){
      nodes.x <- which(x == i)
      nodes.y <- which(y == j)
      
      COST[i,j] <- intersect(nodes.y, nodes.x) %>%
        length()
    }
  }
  
  COST <- 1 / (COST + 1)
  
  perm <- RcppHungarian::HungarianSolver(COST)$pairs[,2]
  
  perm.x <- rep(0, length(x))
  for (i in 1:kmax)
    perm.x[x == i] <- perm[i]
  
  return(perm.x)
}