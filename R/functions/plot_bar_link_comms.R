plot.bar.link.comms <- function(k, net, merde, cluster, labels, regions, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  library(magrittr)
  source("functions/df_to_adj.R")
  source("functions/assign_commship.R")
  source('functions/linkcomm_parameters.R')
  
  tr <- ape::as.phylo(merde)
  tr <- ggtree::fortify(tr)
  tr <- subset(tr, isTip)
  tr <- tr[order(tr$angle, decreasing = T),]
  tr.labels <- tr$label
  nt <- max(net$target)
  ns <- max(net$source)
  min.n <- min(c(nt, ns))
  net <- assign.commship(net, cluster, k,  0, kactive=T)
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  net.mat <- with(net[net$source <= min.n & net$target <= min.n,], data.frame(source=source, target=target, weight=commship)) %>%
    df.to.adj()
  commships <- unique(net.mat[net.mat > 0]) %>%
    sort()
  nc <- commships %>% length()
  t.count <- matrix(0, nrow=min.n, ncol=nc)
  s.count <- matrix(0, nrow=min.n, ncol=nc)
  for (t in 1:min.n){
    for (s in 1:min.n){
      if (net.mat[s,t] >  0){
        t.count[t, which(commships == net.mat[s,t])] <-  t.count[t, which(commships == net.mat[s,t])] + 1
        s.count[t, which(commships == net.mat[t,s])] <-  s.count[t, which(commships == net.mat[t,s])] + 1
      }
    }
  }
  
  t.count <- t.count[match(tr.labels, labels[1:min.n]),]
  s.count <- s.count[match(tr.labels, labels[1:min.n]),]
  data <- data.frame()
  for (e in 1:nc){
    data <- data %>% rbind(data.frame(AREA=tr.labels, COUNT=t.count[,e], TYPE="target", COMMSHIP=commships[e]))
    data <- data %>% rbind(data.frame(AREA=tr.labels, COUNT=s.count[,e], TYPE="source", COMMSHIP=commships[e]))
  }
  
  order.x <- order(regions$COLOR)
  regions$AREA <- regions$AREA %>% tolower()
  data$AREA <- factor(data$AREA, levels = regions$AREA[order.x])
  data$COMMSHIP <- factor(data$COMMSHIP, levels = commships)

  
  p <- ggplot2::ggplot(data, ggplot2::aes(AREA, COUNT, fill=COMMSHIP))+
    ggplot2::facet_grid(TYPE~., scales = "free")+
    ggplot2::geom_bar(color=rgb(0.5,0.5,0.5,0.4), alpha=0.7, stat="identity")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 13, color=regions$COLOR[order.x]))
  
  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 20, height = 10, units = 'in', res = 200)
  print(p)
  dev.off()
}