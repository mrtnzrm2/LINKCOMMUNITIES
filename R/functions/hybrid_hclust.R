library(magrittr)

nodes <- function(merge, members){
  nodes <- list(merge=merge, members=members)
  attr(nodes, "class") <- 'nodes'
  nodes
}

hybrid_hclust <- function(dismat, N, th){

  merge <- matrix(0, nrow = N-1, ncol = 2)
  height <- matrix(0, nrow = N-1, ncol = 1)
  merging.list <- vector(mode = "list", length = N)
  
  for (i in 1:length(merging.list)){
    merging.list[[i]] <- nodes(-i,i)
    
  }
  
  for (j in 1:(N-1)){
    print(j)
    dis <- Inf
    list.length <- N - j + 1
    combos <- combn(list.length, 2)
    
    for(i in 1:ncol(combos)){
      if (j <= th){
        tmp <- min(dismat[merging.list[[combos[1,i]]]$members, merging.list[[combos[2,i]]]$members], na.rm = T)
        if (tmp < dis){
          idx1 <- merging.list[[combos[1,i]]]$merge
          idx2 <- merging.list[[combos[2,i]]]$merge
          cidx1 <- combos[1,i]
          cidx2 <- combos[2,i]
          dis <- tmp
        }
      } else{
        # tmp <- max(dismat[merging.list[[combos[1,i]]]$members, merging.list[[combos[2,i]]]$members], na.rm = T)
        tmp <- mean(dismat[merging.list[[combos[1,i]]]$members, merging.list[[combos[2,i]]]$members], na.rm = T)
        if (tmp < dis){
          idx1 <- merging.list[[combos[1,i]]]$merge
          idx2 <- merging.list[[combos[2,i]]]$merge
          cidx1 <- combos[1,i]
          cidx2 <- combos[2,i]
          dis <- tmp
        }
      }
    }
    
    merge[j,1] <- idx1
    merge[j,2] <- idx2
    
    merging.list[[cidx1]]$members <- c(merging.list[[cidx1]]$members, merging.list[[cidx2]]$members)
    merging.list <- merging.list[-cidx2]
    merging.list[[cidx1]]$merge <- j
    
    height[j] <- dis/2
    
  }
  x <- list(merge=merge, height=height)
  return(x)
}


