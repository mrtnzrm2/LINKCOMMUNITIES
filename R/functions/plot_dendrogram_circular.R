plot.dendrogram.circular <- function(net.cluster, h = -1, k = -1, filename='', plt = T, scale = '', alternative = T,
                                    tip.labels.size = 3, legend.plot.size = 50,
                                    strip.text.size = 2, offset.strip.text = 0.6,
                                    offset.strip = 0.1, branch.none = F, line.size=3,
                                    polar.text.size = 10, bar.line.size=0.8,
                                    tree.regions=list(F,c()),
                                    foldername='', subfolder='',
                                    animal='monkey', path=''){
  
  source('functions/gg_color_hue.R')
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  max.height <- max(net.cluster$height)
  
  if (h>0){
    cut.tree <- cutree(net.cluster, h=h)
  }
  if (k>0){
    cut.tree <- cutree(net.cluster, k=k)
  }
  tr <- ape::as.phylo(net.cluster)
  
  tr <- tidytree::as_tibble(tr)
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
  if (tree.regions[[1]]){
    
    regions.name <- tr$label
    regions.name <- regions.name[!is.na(regions.name)]
    region.df <-tree.regions[[2]] 
    region.df <- region.df[which(region.df$AREA %in% regions.name),]
    region.df <- region.df[match(regions.name, region.df$AREA),]
    n.regions <- length(unique(region.df$REGION))
    tr$regions <- NA
    tr$regions[which(tr$label %in% region.df$AREA)] <- region.df$REGION[which(tr$label %in% region.df$AREA)]
    
    
    tip.colors <- region.df[which(tr$label %in% region.df$AREA), c('REGION', 'COLOR')] 
    tip.colors <- tip.colors[!duplicated(tip.colors$REGION),]
    tip.colors <- tip.colors$COLOR[order(tip.colors$REGION)]

  }
  
  tr <- tidytree::as.treedata(tr)
  if (alternative){
    if (branch.none){
      tree.plot <- ggtree::ggtree(tr,layout='circular', branch.length = 'none', aes(color = as.factor(cluster)), size=line.size)
    } else{
      tree.plot <-  ggtree::ggtree(tr,layout='circular', ggplot2::aes(color = as.factor(cluster)), size=line.size)
    } 
    tree.plot <- tree.plot +  ggtree::geom_tiplab(offset = offset.strip.text, cex=tip.labels.size)
  } else {
    if (branch.none){
      tree.plot <-  ggtree::ggtree(tr,layout='circular', branch.length = 'none', color='blue', size=line.size)
    } else{
      tree.plot <-  ggtree::ggtree(tr,layout='circular', color='blue', size=line.size)
    }
    
    tree.plot <- tree.plot +  ggtree::geom_tiplab(offset = offset.strip.text, cex=tip.labels.size, ggplot2::aes(color = as.factor(cluster)))
  }
  tree.plot <- tree.plot +  ggplot2::scale_color_manual(values = c(gg.color.hue(n-1),'gray'))
  
  if (tree.regions[[1]])
    tree.plot <- tree.plot + 
    ggtree::geom_tippoint(ggplot2::aes(fill=regions), size=7, shape=22)+
    ggplot2::scale_fill_manual(values = tip.colors)
  # scale_fill_brewer(length(n.regions), palette = 'Dark2')
  
  if (plt){
    png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 12, height = 10, units = 'in', res = 200)
    print(tree.plot)
    dev.off()
  }
  else{
    return(tree.plot)
  }
}