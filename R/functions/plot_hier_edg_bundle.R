plot.hier.edg.bundle <- function(K, net, net.cluster, regions, labels, membership, path="", foldername="", subfolder="", filename=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  source('functions/assign_commship.R')
  source('functions/linkcomm_parameters.R')
  source("functions/gg_color_hue.R")
  
  # print("** Using hclust as passed")
  # source('functions/assign_commship.R')
  # net <- assign.commship(net, hcluster, k, 0)
  
  print("*** Warning: You are using a reference hclust to munkres the current hclust")
  source('functions/assign_commship_reference.R')
  net <- assign.commship.reference(net, hcluster, k)
  
  est.para <- linkcomm.parameters(net)
  net$commship[which(net$commship %in% est.para$commship[which(est.para$Dc <= 0)])] <- -1
  net$commship[which(net$commship %in% est.para$commship[is.na(est.para$Dc)])] <- -1
  net <- net[net$commship > 0,]
  net$slabel <- labels[net$source]
  net$tlabel <- labels[net$target]
  nc <- membership %>% unique() %>% length()
  hier <- data.frame(from="origin", to=paste("group", 1:nc, sep = "_"))
  for (cs in 1:nc){
    hier <- hier %>% rbind(data.frame(from=paste("group", cs, sep = "_"), to=labels[membership == cs]))

  }
  
  vertices <- data.frame(name=unique(c(hier$from, hier$to)))
  hier.graph <- igraph::graph_from_data_frame(hier, vertices=vertices)
  hier.vertices <- igraph::V(hier.graph)$name
  hier.labels = hier.vertices[hier.vertices %in% labels]
  
  from <- match(net$slabel, vertices$name)
  to <- match(net$tlabel, vertices$name)

  p <- ggraph::ggraph(hier.graph, layout = 'dendrogram', circular = TRUE) +
    ggraph::geom_conn_bundle(data = ggraph::get_con(from = from, to = to, col=as.factor(net$commship)), ggplot2::aes(color=col), alpha=0.2, width=0.5, tension = 0.9) +
    ggraph::scale_edge_color_manual(values = gg.color.hue(K))+
    ggraph::geom_node_point(ggplot2::aes(filter = leaf, x = x*1.05, y=y*1.05)) +
    ggraph::geom_node_text(ggplot2::aes(x = x*1.1, y=y*1.1, filter = leaf, label=c(NA, rep(NA, nc), hier.labels), angle = ifelse(ggraph::node_angle(x,y) > 90 & ggraph::node_angle(x,y) < 270, ggraph::node_angle(x,y) + 180, ggraph::node_angle(x,y) ), hjust='outward'),
                           size=5, alpha=1, fontface = "bold", color=regions$COLOR[match(hier.labels, regions$AREA)]) +
    ggplot2::theme_void()+
    ggplot2::theme(legend.position = "None")

  png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 12.5, height = 12, units = 'in', res = 200)
  print(p)
  dev.off()
}