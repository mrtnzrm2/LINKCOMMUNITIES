arrange.membership.manual <- function(id){
  # 3 -> 1 1 -> 3 5 -> 6 6 -> 5
  for (i in 1:length(id)){
    if (id[i] == 1)
      id[i] <- 3
    else if (id[i] == 3)
      id[i] <- 1
    else if (id[i] == 5)
      id[i] <- 6
    else if (id[i] == 6)
      id[i] <- 5
  }
  
  return(id)
}

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

plot.heatmap <- function(net, labels){
  print("Warning: Be careful choosing the right partition")
  node.membership <- read.csv("../WSBM/CD/CSV/labels/merged/tracto2016/zz_model/4_r_6_3.csv", header = F) %>% as.matrix()
  node.membership <- arrange.membership.manual(node.membership)
  bourders <- get.bourder.labels(node.membership, labels)
  order.nodes <- order(node.membership)
  net$source.labels <- labels[net$source]
  net$target.labels <- labels[net$target]
  tr.labels <- labels[order.nodes]
  A <- with(net, data.frame(source=source.labels, target=target.labels, weight=log10(fln)+7))
  A$source <- factor(A$source , levels = rev(tr.labels))
  A$target <- factor(A$target, levels = tr.labels)
  p <- ggplot2::ggplot(A, ggplot2::aes(target, source, fill=weight))+
    ggplot2::geom_raster(hjust = 0, vjust = 1)+
    viridis::scale_fill_viridis(option = "A")+
    ggplot2::theme_classic()+
    ggplot2::geom_vline(xintercept = bourders, color="#3DED97", size=1)+
    ggplot2::geom_hline(yintercept = bourders, color="#3DED97", size=1)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size = 8),
                   axis.text.y = ggplot2::element_text(size = 8))
  png("../CSV/merged/imputation/tracto2016/zz_model/fln_4_r_6_3.png", width = 12, height = 10, res = 200, units = "in")
  print(p)
  dev.off()
}