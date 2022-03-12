format.lincomm <- function(net, hcluster, k){
  source("functions/assign_commship_reference.R")
  net <- assign.commship.reference(net, hcluster, k)
  net$commship[net$commship == 3] <- 1
  net$commship[net$commship == 4] <- 3
  return(net)
}