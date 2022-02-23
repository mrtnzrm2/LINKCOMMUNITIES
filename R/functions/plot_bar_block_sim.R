get.block.commship <- function(A, R, wsbm){
  block.commship <- matrix(0, R, R)
  
  for (i in 1:R){
    if (i == 1){
      x.i <- 1
      x.f <- wsbm[1]
    } else{
      x.i <- sum(wsbm[1:(i-1)]) + 1
      x.f <- sum(wsbm[1:i])
    }
    for (j in 1:R){
      if (j == 1){
        y.i <- 1
        y.f <- wsbm[1]
      } else{
        y.i <- sum(wsbm[1:(j-1)]) + 1
        y.f <- sum(wsbm[1:j])
      }
      
      commships <- A[x.i:x.f,y.i:y.f] %>% as.array() %>% table()
      if ("0" %in% names(commships))
        commships <- commships[names(commships) != "0"]
      nc <- commships %>% names()
      
      block.commship[i,j] <- nc[commships == max(commships)] %>% as.numeric()
    }
  }
  return(block.commship)
}

filter.block.commship <- function(k, K, data){
  
  data <- data[data$source == k | data$target == k,]
  df <- data.frame(u=rep(k, K), v=1:K, weight=rep(0, K), count=rep(0, K))
  
  for (i in 1:nrow(data)){
    if (data$source[i] == k){
      v <- data$target[i]
    } else{
      v <- data$source[i]
    }
    df$weight[df$v == v] <- data$sim[i] + df$weight[df$v == v] 
    df$count[df$v == v] <- df$count[df$v == v]  + 1
  }
  
  for (i in 1:K){
    if (df$count[i] != 0){
      df$weight[i] <- df$weight[i]/df$count[i]
    } else
      df$weight[i] <-  0
  }
  
  return(df)
}

plot.bar.block.sim <- function(R, K, net, net.cluster, path="", foldername="", subfolder="", filename=""){
  
  source("functions/assign_commship.R")
  source("functions/linkcomm_parameters.R")
  source("functions/df_to_adj.R")
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  print("Warning: Be careful choosing the right partition")
  wsbm.partition <- read.csv("../WSBM/CD/CSV/labels/merged/tracto2016/zz_model/4_r_6_3.csv", header = F) %>% as.matrix()
  order.wsbm <- order(wsbm.partition)
  table.wsbm <- table(wsbm.partition)
  
  net <- assign.commship(net, net.cluster, K, 0)
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  
  A <- with(net, data.frame(source=source, target=target, weight=commship)) %>% df.to.adj()
  A <- A[order.wsbm, order.wsbm]
  A <- get.block.commship(A, R, table.wsbm)
  
  print("Warning: Load correct node and sim blocks")
  nodes.list <- read.csv("../CSV/merged/similarity/tracto2016/zz_model/nodes_block_sim_4_r_6_3.csv", header = F) %>% as.matrix()
  sim.list <- read.csv("../CSV/merged/similarity/tracto2016/zz_model/sim_block_sim_4_r_6_3.csv", header = F) %>% as.matrix() %>% as.numeric()
  data <- data.frame(u.x=nodes.list[,1], u.y=nodes.list[,2], v.x=nodes.list[,3], v.y=nodes.list[,4], sim=sim.list)
  data$label.u <- 0
  data$label.v <- 0
  data$beer < ""
  data$source <- 0
  data$target <- 0
  
  for (e in 1:nrow(data)){
    data$label.u[e] <- R*(data$u.x[e])+ data$u.y[e] + 1
    data$label.v[e] <- R*(data$v.x[e])+ data$v.y[e] + 1
    data$source[e] <- A[data$u.x[e]+1, data$u.y[e]+1]
    data$target[e] <- A[data$v.x[e]+1, data$v.y[e]+1]
    if (data$label.u[e] > data$label.v[e]){
      data$beer[e] <- paste0(data$label.v[e],data$label.u[e]) 
    }
    else{
      data$beer[e] <- paste0(data$label.u[e],data$label.v[e]) 
    }
  }
  
  data <- data[!duplicated(data$beer),]
  data.block <- data.frame()
  
  K.in.A <- pracma::Reshape(A, R*R, 1) %>% unique() %>% sort()
  for (k in K.in.A)
    data.block <- data.block %>% rbind(filter.block.commship(k, K, data))
  
  p <- ggplot2::ggplot(data.block, ggplot2::aes(v, weight))+
    ggplot2::facet_wrap(~u)+
    ggplot2::geom_bar(stat = "identity")
  
  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 12, height = 10, units = 'in', res = 200)
  print(p)
  dev.off()
}