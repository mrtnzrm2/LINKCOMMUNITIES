plot.dendrogram.matrix <- function(net, net.merde, net.cluster, labels, regions, nodes,
                                  h.merde = -1, k.merde = -1, k=-1, height=0, filename='', plt = T, 
                                  alternative = T, tip.labels.size = 3, 
                                  branch.none = F, line.size=3, foldername='', 
                                  subfolder='', path=''){
  
  source('functions/df_to_adj.R')
  source('functions/extract_k_partition.R')
  source('functions/with_potholes.R')
  source('functions/gg_color_hue.R')
  source('functions/adj_to_df.R')
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.height <- max(net.merde$height)
  net <- net[,c('source', 'target', 'weight', 'id')]
  A <- df.to.adj(net[1:3])
  A[A == 0] <- NA
  A <- log(A) + 7
  A <- A %>% as.data.frame()
  colnames(A) <- labels
  rownames(A) <- labels
  
  if (h.merde>0){
    cut.tree <- cutree(net.merde, h=h.merde)
  }
  if (k.merde>0){
    cut.tree <- cutree(net.merde, k=k.merde)
  }
  
  tr <- ape::as.phylo(net.merde)
  
  trd <- ggtree::fortify(tr)
  trd <- subset(trd, isTip)
  trd <- trd[order(trd$angle, decreasing = T),]
  
  historia <- extract.k.partition(net, net.cluster, nodes, k, height)
  historia <- historia[match(trd$label, labels[1:nodes])]
  historia <- with.potholes(historia)
  
  u.hist <- historia[[1]]
  tab.t <- historia[[2]]
  tab2.t <- tab.t
  
  for (e in 2:length(tab.t)){
    tab.t[e] <- tab.t[e] + tab.t[e-1]
  }
  tab.t <- tab.t + 0.5
  tab.op <- c(0, tab.t)
  for (e in 1:length(tab2.t)){
    tab2.t[e] <- tab.t[e] - (tab.op[e+1] - tab.op[e])/2
  }
  
  tab2.t <- as.integer(tab2.t)
  
  trd.names <- as.character(trd$label)
  tr <- tibble::as_tibble(tr)
  
  n <- length(unique(cut.tree))
  anc <- matrix(0, nrow = n, ncol = 1)
  anc.tips <- matrix(0, nrow = n, ncol = 2)
  ss <- cut.tree[net.cluster$order]
  sort.ss <- sort(unique(cut.tree))
  tg <- tibble::tibble('node' = c(0), 'cluster'= c(0))
  
  if (alternative){
    for (i in 1:n){
      sp <- tr$node[which(cut.tree == sort.ss[i])]
      anc.tips[i,1] <- sp[1]
      anc.tips[i,2] <- sp[length(sp)]
      anc[i] <- ggtree::MRCA(tr, sp)$node
    }
    
    tg <- tibble::tibble( 'node' = c(0), 'cluster' = c(0))
    for (i in 1:n){
      offsp <- tidytree::offspring(tr, anc[i])$node
      tg <- dplyr::full_join(tg, tibble::tibble('node' = offsp, 'cluster' = rep(sort.ss[i],length(offsp))), by = c('node', 'cluster'))
    }
  } else {
    for (i in 1:n){
      sp <- tr$node[which(cut.tree == sort.ss[i])]
      tg <- dplyr::full_join(tg, tibble::tibble('node' = sp, 'cluster' = sort.ss[i]), by = c('node', 'cluster'))
      anc.tips[i,1] <- sp[1]
      anc.tips[i,2] <- sp[length(sp)]
    }
  }
  tg <- tg[2:nrow(tg),]
  tr <- dplyr::full_join(tr, tg, by = 'node')
  tr$cluster[is.na(tr$cluster)] <- 'NA'
  n <- length(unique(tr$cluster))
  
  region.df <- regions[which(regions$AREA %in% trd.names),]
  region.df <- region.df[match(trd.names, region.df$AREA),]
  tr$regions <- NA
  tr$regions[which(tr$label %in% region.df$AREA)] <- region.df$REGION[which(tr$label %in% region.df$AREA)]
  color.regions <- regions[order(regions$REGION),]
  
  tr <- tidytree::as.treedata(tr)
  if (alternative){
    if (branch.none){
      tree.plot <- ggtree::ggtree(tr,layout='rectangular', branch.length = 'none', ggplot2::aes(color = as.factor(cluster)), size=line.size)
    } else{
      tree.plot <- ggtree::ggtree(tr,layout='rectangular', ggplot2::aes(color = as.factor(cluster)), size=line.size)
    } 
    tree.plot <- tree.plot #+ geom_tiplab(offset = 0.03, cex=tip.labels.size)
  } else {
    if (branch.none){
      tree.plot <- ggtree::ggtree(tr,layout='rectangular', branch.length = 'none', color='blue', size=line.size)
    } else{
      tree.plot <- ggtree::ggtree(tr,layout='rectangular', color='blue', size=line.size)
    }
    
    tree.plot <- tree.plot #+ geom_tiplab(cex=tip.labels.size, aes(color = as.factor(cluster)))
  }
  
  tree.plot <- tree.plot + 
    ggplot2::scale_color_manual(values = c(gg.color.hue(n-1), 'grey'))+
    ggplot2::guides(color=ggplot2::guide_legend(title="Link community"))
  
  leg.t <- cowplot::get_legend(tree.plot)
  
  A <- A[trd.names, trd.names] %>%
    adj.to.df()
  A <- A[!is.na(A$weight),]
  
  A <- ggplot2::ggplot(A, ggplot2::aes(target, source))+
    ggplot2::geom_tile(ggplot2::aes(fill=weight))+
    ggplot2::scale_x_continuous(expand = c(0, 0))+
    ggplot2::scale_y_continuous(trans = 'reverse',expand = c(0, 0))+
    viridis::scale_fill_viridis(direction = 1, option = 'A')
  
  A <- A + ggplot2::geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(0.5,0.5)), ggplot2::aes(x,y), color='green', size=1)+
    ggplot2::geom_line(data=data.frame(x=c(0.5,0.5), y=c(0.5,tab.t[1])), ggplot2::aes(x,y), color='green', size=1)+
    ggplot2::geom_line(data=data.frame(x=c(tab.t[1],tab.t[1]), y=c(0.5,tab.t[1])), ggplot2::aes(x,y), color='green', size=1)+
    ggplot2::geom_line(data=data.frame(x=c(0.5,tab.t[1]), y=c(tab.t[1],tab.t[1])), ggplot2::aes(x,y), color='green', size=1)
  
  xp <- tab.t[1]
  yp <- xp
  
  for (e in 2:(length(tab.t)-1)){
    A  <- A + ggplot2::geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(yp,yp)), ggplot2::aes(x,y), color='green', size=1)+
      ggplot2::geom_line(data=data.frame(x=c(xp,xp), y=c(yp,tab.t[e])), ggplot2::aes(x,y), color='green', size=1)+
      ggplot2::geom_line(data=data.frame(x=c(tab.t[e],tab.t[e]), y=c(yp,tab.t[e])), ggplot2::aes(x,y), color='green', size=1)+
      ggplot2::geom_line(data=data.frame(x=c(xp,tab.t[e]), y=c(tab.t[e],tab.t[e])), ggplot2::aes(x,y), color='green', size=1)
    
    xp <- tab.t[e]
    yp <- xp
  }
  
  
  leg.a <- cowplot::get_legend(A)
  
  A <- A + ggplot2::guides(fill=ggplot2::guide_legend(title = 'w'))+
    ggplot2::theme_void()
  
  region.df$position <- 1
  region.df$NODE <- 1:nodes
  
  gg_leg_y <- ggplot2::ggplot(region.df, ggplot2::aes(y=NODE, x=position, fill=REGION))+
    ggplot2::geom_tile()+
    ggplot2::scale_x_continuous(expand = c(0, 0))+
    ggplot2::scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
    # scale_fill_manual(values = unique(color.regions$COLOR))+
    ggplot2::geom_text(label=trd.names, color='white')+
    ggplot2::guides(fill=ggplot2::guide_legend(title = 'Area'))+
    ggplot2::theme_void()
  
  leg.gg_y <- cowplot::get_legend(gg_leg_y)
  
  gg_leg_x <- ggplot2::ggplot(region.df, ggplot2::aes(x=NODE, y=position, fill=REGION))+
    ggplot2::geom_tile()+
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0, 0))+
    # scale_fill_manual(values = unique(color.regions$COLOR))+
    ggplot2::geom_text(label=trd.names, color='white', angle=90)+
    ggplot2::guides(fill=ggplot2::guide_legend(title = 'Area'))+
    ggplot2:: theme_void()
  
  leg.gg_x <- cowplot::get_legend(gg_leg_x)
  
  # gg_leg <- ggplot(region.df, aes(y=NODE, x=position, fill=REGION))+
  #   geom_tile()+
  #   scale_x_continuous(expand = c(0, 0))+
  #   # geom_hline(yintercept = tab.t[1:(length(tab.t)-1)], color='pink', size=1)+
  #   scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
  #   # scale_fill_manual(values = unique(color.regions$COLOR))+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   theme_void()
  # 
  # leg.gg <- get_legend(gg_leg)
  
  legend <- cowplot::plot_grid(leg.t, leg.a, leg.gg_y, nrow = 3, rel_heights = c(2,1,5), axis = 't')
  
  # g <- plot_grid(tree.plot + theme(legend.position = 'none'), 
  #                A + theme(legend.position = 'none',
  #                          panel.grid = element_blank(),
  #                          panel.border = element_blank()), 
  #                gg_leg + theme(legend.position = 'none',
  #                               panel.grid = element_blank(),
  #                               panel.border = element_blank()), 
  #                legend, ncol=4, rel_widths=c(5,15,1,2), 
  #                align = 'h' ,
  #                scale=0.98)
  
  g <- cowplot::plot_grid(NULL, gg_leg_x + ggplot2::theme(legend.position = 'none',
                                        panel.grid = ggplot2::element_blank(),
                                        panel.border = ggplot2::element_blank()), 
                 NULL, NULL,
                 tree.plot + ggplot2::theme(legend.position = 'none'),
                 A + ggplot2::theme(legend.position = 'none',
                           panel.grid = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank()),
                 gg_leg_y + ggplot2::theme(legend.position = 'none',
                                  panel.grid = ggplot2::element_blank(),
                                  panel.border = ggplot2::element_blank()), 
                 legend, 
                 ncol=4, nrow = 2, rel_widths=c(5,15,1,2), rel_heights = c(1,15),
                 align = 'hv' ,
                 scale=0.98)
  
  if (plt){
    png( sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 20, height = 13, units = 'in', res = 200)
    print(g)
    dev.off()
  }
  else{
    print(tree.plot)
  }
}