cascade.arc.en.ciel <- function(net.merde, nodes, labels, regions, cascade, plt=F, subfolder='',
                                filename='', path='', foldername='', animal=''){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  tr <- ape::as.phylo(net.merde)
  
  trd <- ggtree::fortify(tr)
  trd <- subset(trd, isTip)
  trd <- trd[order(trd$angle, decreasing = T),]
  trd.names <- trd$label
  order.trd <- match(trd.names, labels[1:nodes])
  
  regions <- regions[which(regions$AREA %in% labels[1:nodes]),]
  regions <- regions[match(labels[1:nodes], regions$AREA),]
  
  cascade <- cascade[,order.trd]
  df.cascade <- adj.to.df(cascade)
  df.cascade$weight <- as.factor(df.cascade$weight)
  
  regions <- regions[order.trd,]
  regions$NODE <- 1:nodes
  regions$position <- 1
  
  gg_leg <- ggplot2::ggplot(regions, ggplot2::aes(y=position, x=NODE, fill=REGION))+
    ggplot2::geom_tile()+
    ggplot2::scale_x_continuous(expand = c(0, 0))+
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    # scale_fill_manual(values = unique(color.regions$COLOR))+
    ggplot2::guides(fill=ggplot2::guide_legend(title = 'Area'))+
    ggplot2::annotate(geom = "text", x=regions$NODE, y=regions$position,
             label = regions$AREA, color = "white",
             angle = 90, fontface=2, size=2)+ 
    ggplot2::theme_void()
  
  leg.gg <- cowplot::get_legend(gg_leg)
  
  pp <- ggplot2::ggplot(df.cascade, ggplot2::aes(target, source, fill=weight))+
    ggplot2::geom_tile()+
    ggplot2::scale_fill_discrete(na.value='black')+
    ggplot2::scale_x_continuous(expand = c(0, 0))+
    ggplot2::scale_y_continuous(trans = 'reverse', expand = c(0, 0),
                       breaks = 1:nodes, labels = 1:nodes)+
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size=3),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x=ggplot2::element_blank(),
      axis.title.y=ggplot2::element_blank())
  
  pat.1 <- cowplot::plot_grid(gg_leg + ggplot2::theme(legend.position = 'none'),
                     pp + ggplot2::theme(legend.position = 'none'),
                     nrow = 2, rel_heights =c(1,10), 
                     align = 'v' ,
                     scale=0.98)
  
  pate <- cowplot::plot_grid(pat.1, leg.gg,
                    ncol = 2, rel_widths =c(12, 1), 
                    align = 'h' ,
                    scale=0.98, axis = 'l')
  
  if (plt){
    png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 10, height = 7, units = 'in', res = 200)
    print(pate)
    dev.off()
  }
  else{
    print(pate)
    # plot_gg(pate, width = 20, height = 13, raytrace = F, preview = TRUE)
  }
}