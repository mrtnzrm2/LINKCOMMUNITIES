plot.densities <- function(net, labels, memberships){
  source("functions/get_bourder_labels.R")

  bourders <- get.bourder.labels(memberships, labels)
  order.nodes <- order(memberships)
  memberships <- memberships[order.nodes]
  A <- with(net, data.frame(source=source, target=target, weight=fln))
  source("functions/df_to_adj.R")
  A <- A %>% df.to.adj()
  A <- A[order.nodes, order.nodes]
  membs <- memberships %>% unique() %>% sort()
  nmbs <- membs %>% length() 
  densities <- matrix(0, nrow = nmbs, ncol = nmbs)
  for (i in membs){
    for (j in membs){
      if (i == j){
        den.i <- which(memberships == i)
        n.i <- den.i %>% length()
        a <- A[den.i,den.i]
        densities[i,j] <- sum(a > 0)/(n.i*(n.i - 1))
      } else{
        den.i <- which(memberships == i)
        den.j <- which(memberships == j)
        n.i <- den.i %>% length()
        n.j <- den.j %>% length()
        a <- A[den.i,den.j]
        densities[i,j] <- sum(a > 0)/(n.i*n.j)
      }
    }
  }
  print(densities)
  A <- matrix(0, nrow = 107, ncol = 107)
  source("functions/adj_to_df.R")
  A <- adj.to.df(A)
  A$source.membership <- memberships[A$source]
  A$target.membership <- memberships[A$target]
  A$density <- 0
  for (e in 1:nrow(A)){
    A$density[e] <- densities[A$target.membership[e], A$source.membership[e]]
  }
  tr.labels <- labels[order.nodes]
  A$source.label <- tr.labels[A$source]
  A$target.label <- tr.labels[A$target]
  A$source.label <- factor(A$source.label, levels = tr.labels)
  A$target.label <- factor(A$target.label, levels = rev(tr.labels))
  p <- ggplot2::ggplot(A, ggplot2::aes(source.label, target.label, fill=density))+
    ggplot2::geom_raster(hjust = 0, vjust = 1)+
    viridis::scale_fill_viridis(option = "C", direction = 1)+
    ggplot2::geom_hline(yintercept = bourders, size=0.5)+
    ggplot2::geom_vline(xintercept = bourders, size=0.5)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size = 4),
                   axis.text.y = ggplot2::element_text(size=4))
  
  png("../CSV/merged/imputation/tracto2016/zz_model/densities_4_r_6_3.png", width = 6, height = 5, res = 200, units = "in")
  print(p)
  dev.off()
  
}
