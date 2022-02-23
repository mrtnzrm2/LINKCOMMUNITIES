linkcomm.parameters <- function(net){
  net <- net[order(net$id),]
  
  num.comm <- length(unique(net$commship))
  comm <- unique(net$commship)
  
  edense <- rep(0, num.comm)
  clustering <- rep(0, num.comm)
  Dc <- rep(NA, num.comm)
  n <- rep(0, num.comm)
  
  for (i in 1:num.comm){
    ss <- net[which(net$commship == comm[i]), c('source','target','weight', 'id')]
    n[i] <- length(unique(c(ss$source, ss$target)))
    eds <- length(ss$weight)
    edense[i] <- eds/(n[i]*(n[i]-1))
    clustering[i] <- sum(log10(ss$weight) + 7)/sum(log10(net$weight) + 7)
    
    leaves <- nrow(ss)
    Dc[i] <- (leaves - n[i] + 1)/(n[i]-1)^2
  }
  est.data <- data.frame('EdgeDensity'= edense, 
                         'W.clustering' = clustering,
                         'commship' = comm, 
                         'Nodes'= n,
                         'Dc' = Dc)
  
  est.data <- est.data[order(est.data$Nodes, decreasing=T),]
  return(est.data)
}
