create.id.matrix <- function(A, N){
  print("*** Matrix has to be squared for the similarity competitive function to work")
  id <- matrix(0, ncol = N, nrow = N)
  count <- 1
  for (i in 1:N){
    for (j in 1:N){
      if (A[i,j] > 0){
        id[i,j] <- count
        count <- count + 1
        }
    }
  }
  return(id)
}

similarity.competitive <- function(net, leaves, self.loop = F, mode=''){
  source('functions/jaccard_p_fast.R')
  source("functions/aki.R")
  source("functions/aik.R")
  source("functions/df_to_adj.R")
  
  aki.matrix <- aki(net, self.loop=self.loop, mode=mode)
  aik.matrix <- aik(net, self.loop=self.loop, mode=mode)
  
  net <- net %>% df.to.adj()
  N <- nrow(net)
  id <- create.id.matrix(net, N)
  net.sim <- matrix(0, nrow = leaves, ncol = leaves)
  
  for (i in 1:N){
    for (j in 1:N){
      id.1 <- id[i,j]
      if (id.1 <= 0)
        next
      
      for (k in 1:N){
        id.2 <- id[i,k]
        if (id.2 <= 0)
          next
        net.sim[id.1, id.2] <- jaccard.p.fast(aki.matrix[j,], aki.matrix[k,])
      }
      for (k in 1:N){
        id.2 <- id[k,j]
        if (id.2 <= 0)
          next
        net.sim[id.1, id.2] <- jaccard.p.fast(aik.matrix[i,], aik.matrix[k,])
      }
    }
  }
  net.sim <- net.sim + t(net.sim)
  net.sim[net.sim == 0] <- NA
  return(net.sim)
}