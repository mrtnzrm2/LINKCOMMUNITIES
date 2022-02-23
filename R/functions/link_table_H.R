link.table.H <- function(net, net.cluster, labels, k){
  library(magrittr)
  source('functions/assign_commship.R')
  source('functions/linkcomm_parameters.R')
  
  net <- assign.commship(net, net.cluster, k, 0)
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  
  ncol.areas <- max(net$target)
  nrow.areas <- max(net$source)
  
  plink <- data.frame(AREA=c(), DIC=c(), COMMSHIP=c(), P=c())
  for (i in 1:ncol.areas){
    v.commships <- net[net$target == i, c("commship", "weight")]
    H <- sum(v.commships$weight)
    commships <- v.commships$commship %>% unique() %>% sort()
    nc <- commships %>% length()
    for (j in 1:nc){
      cwsum <- v.commships$weight[v.commships$commship == commships[j]] %>% sum()
      plink <- plink %>% rbind(data.frame(AREA=labels[i], DIC='V', COMMSHIP=commships[j], P=cwsum/H)) 
    }
  }
  for (i in 1:nrow.areas){
    h.commships <- net[net$source == i, c("commship", "weight")]
    H <- sum(h.commships$weight)
    commships <- h.commships$commship %>% unique() %>% sort()
    nc <- commships %>% length()
    for (j in 1:nc){
      cwsum <- h.commships$weight[h.commships$commship == commships[j]] %>% sum()
      plink <- plink %>% rbind(data.frame(AREA=labels[i], DIC='H', COMMSHIP=commships[j], P=cwsum/H)) 
    }
  }
  return(plink)
}