get.bourder.labels <- function(memberships, labels){
  library(magrittr)
  get.order <- order(memberships)
  memberships <- memberships[get.order]
  labels <- labels[get.order]
  nm <- memberships %>% unique() %>% length()
  tab.mem <- memberships %>% table()
  bourders <- rep("", nm-1)
  for (i in 1:nm){
    bourders[i] <- labels[sum(tab.mem[1:i])]
  }

  return(bourders[1:(length(bourders)-1)])
}

plot.lincomm.matrix <- function(net, net.cluster, 
                                leaves, nodes, k, best.height, labels, 
                                regions, plt=F,
                                filename='', foldername='', subfolder='',
                                path=''){
  library(magrittr)
  source('functions/extract_k_partition.R')
  source('functions/with_potholes.R')
  source('functions/assign_commship.R')
  source('functions/linkcomm_parameters.R')
  source("functions/gg_color_hue.R")
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  # tr <- ape::as.phylo(net.merde)
  # tr <- ggtree::fortify(tr)
  # tr <- subset(tr, isTip)
  # tr <- tr[order(tr$angle, decreasing = T),]
  
  net <- assign.commship(net, net.cluster, k, best.height)
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  commships <- net$commship %>% unique()
  ncommships <- commships %>% length()
  color.palette <- gg.color.hue(ncommships)
  color.commship <- data.frame()
  for (e in 1:ncommships){
    if (commships[e] > 0){
      color.commship <- color.commship %>% rbind(data.frame(commship=commships[e], color=color.palette[e]))
    }
    else{
      color.commship <- color.commship %>% rbind(data.frame(commship=commships[e], color=rgb(0.5,0.5,0.5,0.3))) 
    }
  }
  color.commship <- color.commship[order(color.commship$commship),]
  
  print("Warning: Be careful choosing the right partition")
  node.membership <- read.csv("../WSBM/CD/CSV/labels/merged/tracto2016/zz_model/4_r_6_3.csv", header = F) %>% as.matrix()
  bourders <- get.bourder.labels(node.membership, labels)
  tr.labels <- labels[order(node.membership)]
  
  net <- with(net, data.frame(source=source, target=target, weight=commship))
  A <- net
  A <- A[which(A$weight != 0),]
  # A$source.label <- labels[A$source]
  # A$target.label <- labels[A$target]
  A$source.label <- tr.labels[A$source]
  A$target.label <- tr.labels[A$target]

  A$source.label <- factor(A$source.label, levels = rev(tr.labels))
  A$target.label <- factor(A$target.label, levels = tr.labels)
  


  q <- ggplot2::ggplot(A, ggplot2::aes(target.label, source.label, fill=as.factor(weight)))+
    ggplot2::geom_tile(color='white')+
    ggplot2::geom_text(label=A$weight, size=1)+
    ggplot2::scale_fill_manual(values = color.commship$color)+
    ggplot2::theme_classic()+
    ggplot2::geom_vline(xintercept = bourders)+
    ggplot2::geom_hline(yintercept = bourders)+
    # ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold", angle = 90),
    #                legend.position = "none",
    #                axis.text.y = ggplot2::element_text(face="bold"))
    ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold", angle = 90, color=regions$COLOR[match(tr.labels, regions$AREA)]),
                   legend.position = "none",
                   axis.text.y = ggplot2::element_text(face="bold", color=regions$COLOR[match(rev(tr.labels), regions$AREA)]))

  if (plt){
    png(sprintf("%s/%s/%s/%s_linkcommunities.png", path, foldername, subfolder, filename), width = 12, height = 10, units = 'in', res = 200)
    print(q)
    dev.off()
  }
  else{
    print(q)
  }
}