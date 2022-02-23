cascade.sil.vous.plaît <- function(net, net.merde, nodes){
  
  library(magrittr)
  source('functions/adj_to_df.R')
  source('functions/df_to_adj.R')
  
  max.sim <- pracma::ceil(-log10(min(net$weight))) + 1
  net$weight <- log10(net$weight) + max.sim
  
  cascade <- matrix(0, ncol = nodes, nrow = nodes)
  cascade[nodes,] <- as.factor(cutree(net.merde, k=nodes))
  
  height <- matrix(1, ncol = nodes, nrow = nodes)
  
  for (e in (nodes-1):1){
    cascade[e,] <- cutree(net.merde, k=e)
    old.cls <- unique(cascade[e+1,])
    new.cls <- unique(cascade[e,])
    
    for (e.o in old.cls){
      for (e.i in new.cls){
        
        new.nodes <- which(cascade[e,] == e.i)
        old.nodes <- which(cascade[e+1,] == e.o)
        
        rep.n <- length(intersect(new.nodes, old.nodes))
        
        if ( rep.n >= 1 && length(new.nodes) > length(old.nodes)){
          height[e,new.nodes] <- height[e+1,new.nodes] + 1
        } else {
          for (n in new.nodes){
            height[e,n] <- max(height[,n])
          }
        }
      }
    }
  }
  
  cascade <- matrix(0, ncol = nodes, nrow = nodes)
  cascade[nodes,] <- 1:nodes
  
  for (e in (nodes-1):1){
    
    cascade[e,] <- cascade[e+1,]
    
    new.coms <- cutree(net.merde, k=e)
    old.coms <- cutree(net.merde, k=e+1)
    
    old.cls <- unique(old.coms)
    new.cls <- unique(new.coms)
    
    for (e.i in new.cls){
      new.nodes <- which(new.coms == e.i)
      n <- 1
      cfirst <- c()
      csecond <- c()
      dog <- unlist(lapply(old.cls, function(x) 
        ifelse(length(intersect(new.nodes, which(old.coms == x))) >= 1,
               T, F)))
      if (sum(dog) == 2){
        cfirst <-  which(old.coms == old.cls[which(dog)[1]])
        csecond <-  which(old.coms == old.cls[which(dog)[2]])
      }
      if (length(cfirst) > 0 && length(csecond) > 0){
        for (n in nodes:e){
          n1 <- sum(height[n, cfirst] > 1)
          n2 <- sum(height[n, csecond] > 1)
          if (n1 > n2 && n1 > 0){
            the.nodes <- cfirst[which(height[n,cfirst] != 1)]
            break
          } else if (n2 > n1 && n2 > 0){
            the.nodes <- csecond[which(height[n,csecond] != 1)]
            break
          } else if (n2 == n1 && n2 > 0){
            the.nodes <- c(cfirst[which(height[n,cfirst] != 1)],
                           csecond[which(height[n,csecond] != 1)])
            break
          }
        }
        cascade[e, new.nodes] <- cascade[n, max(the.nodes)]
      }
    }
  }
  
  # height <- height[,order.trd]
  height <- adj.to.df(height)
  
  # cascade <- cascade[,order.trd]
  df.cascade <- adj.to.df(cascade)
  df.cascade$weight <- as.factor(df.cascade$weight)
  df.cascade$weight[which(height$weight == 1)] <- NA 
  
  cascade <- df.cascade %>% 
    df.to.adj()
  
  return(cascade)
} 