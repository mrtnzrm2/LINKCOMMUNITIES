la.arbre.a.merde <- function(hclust.features, nodes, labels, net, net.cluster, leaves){
  
  source('functions/assign_commship.R')
  source('functions/linkcomm_parameters.R')
  
  k.remark <- hclust.features$K[which(hclust.features$NEC >= 1)]
  loc.max <- hclust.features$height[which(hclust.features$NEC >= 1)]
  loc.dup <- !duplicated(loc.max)
  k.remark <- k.remark[loc.dup]
  loc.max <- loc.max[loc.dup]
  
  me.r <- matrix(0, nrow = nodes - 1, ncol=2)
  height <- rep(0, nodes - 1)
  
  merging.list <- vector(mode = "list", length = nodes)
  for (i in 1:length(merging.list)){
    merging.list[[i]] <- list(-i,i)
  }
  
  ct <- 1
  for (e in 1:length(k.remark)){
    if (mod(e,100) == 0)
      print(e)
    
    k <- k.remark[e]
    best.height <- loc.max[e]
    
    net.reduced <- assign.commship(net[which(net$id <= leaves),], net.cluster, k, best.height, kactive = T)
    est.para <- linkcomm.parameters(net.reduced)
    net.reduced$commship[which(net.reduced$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
    net.reduced$commship[which(net.reduced$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
    net.reduced <- net.reduced[which(net.reduced$commship != -1),]
    
    coms <- sort(unique(net.reduced$commship))
    
    merge.copy <- merging.list
    
    for (com in coms){
      net.com <- net.reduced[which(net.reduced$commship == com),]
      nds <- sort(intersect(unique(net.com$source), unique(net.com$target)))
      if (length(nds) > 1){
        
        dog <- unlist(lapply(1:length(merge.copy), function(x) ifelse(length(intersect(nds, merge.copy[[x]][[2]])) >= 1, T, F)))
        dog.id <- which(dog)
        
        if (length(dog.id) == 2){
          me.r[ct,1] <- merge.copy[[dog.id[1]]][[1]]
          me.r[ct,2] <- merge.copy[[dog.id[2]]][[1]]
          height[ct] <- best.height
          
          merge.copy <- merge.copy[-dog.id]
          
          merge.copy[[length(merge.copy)+1]] <- list(ct, sort(unique(c(merging.list[[dog.id[1]]][[2]], merging.list[[dog.id[2]]][[2]]))))
          
          ct <- ct + 1
        } else if (length(dog.id) > 2){
          
          merge.copy2 <- vector(mode = "list", length = length(dog.id))
          
          for (ii in 1:length(dog.id)){
            merge.copy2[[ii]] <- merge.copy[[dog.id[ii]]]
          }
          
          merge.copy <- merge.copy[-dog.id]
          
          dog2 <- unlist(lapply(1:length(merge.copy2), function(x) ifelse(length(intersect(nds, merge.copy2[[x]][[2]])) >= 1, T, F)))
          dog2.id <- which(dog2)
          
          while(length(merge.copy2) >= 2){
            me.r[ct, 1] <- merge.copy2[[dog2.id[1]]][[1]]
            me.r[ct, 2] <- merge.copy2[[dog2.id[2]]][[1]]
            height[ct] <- best.height
            
            if (length(merge.copy2) > 2){
              merge.copy2[[length(merge.copy2)+1]] <- list(ct, sort(unique(c(merge.copy2[[dog2.id[1]]][[2]], merge.copy2[[dog2.id[2]]][[2]]))))
              merge.copy2 <- merge.copy2[-c(dog2.id[1], dog2.id[2])]
              
              dog2 <- unlist(lapply(1:length(merge.copy2), function(x) ifelse(length(intersect(nds, merge.copy2[[x]][[2]])) >= 1, T, F)))
              dog2.id <- which(dog2)
              
            } else if (length(merge.copy2) == 2){
              
              merge.copy[[length(merge.copy)+1]] <- list(ct, sort(unique(c(merge.copy2[[dog2.id[1]]][[2]], merge.copy2[[dog2.id[2]]][[2]]))))
              merge.copy2 <- merge.copy2[-c(dog2.id[1], dog2.id[2])]
            }
            
            ct <- ct + 1
          }
          
        }
        merging.list <- merge.copy
      }
    }
    
  }
  
  print(merging.list)
  
  merde <- list()
  merde$merge <- me.r
  merde$order <- 1:nodes
  merde$labels <- labels
  merde$height <- height
  class(merde) <- 'hclust'
  
  return(merde)
  
}