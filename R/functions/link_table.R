link.table <- function(net, net.cluster, labels, k){
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
    v.commships <- net$commship[net$target == i] %>% table()
    nv <- names(v.commships)
    for (j in 1:length(v.commships)){
      id.j <- nv[j] %>% as.numeric()
      if (id.j > 0)
        plink <- plink %>% rbind(data.frame(AREA=labels[i], DIC='V', COMMSHIP=id.j, P=v.commships[j]/sum(v.commships))) 
    }
  }
  for (i in 1:nrow.areas){
    h.commships <- net$commship[net$source == i] %>% table()
    nh <- names(h.commships)
    nlc <- 1:length(nh)
    for (j in 1:length(h.commships)){
      id.j <- nh[j] %>% as.numeric()
      if (id.j > 0)
        plink <- plink %>% rbind(data.frame(AREA=labels[i], DIC='H', COMMSHIP=id.j, P=h.commships[j]/sum(h.commships)))
    }
    
  }
  return(plink)
}