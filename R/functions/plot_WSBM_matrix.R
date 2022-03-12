plot.WSBM.matrix <- function(net, labels, membership, path='', subfolder='', filename=''){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s', path, subfolder), showWarnings = FALSE)
  }
  
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/get_bourder_labels.R")
  
  bourders <- get.bourder.labels(membership, labels)
  order.nodes <- order(membership)
  labels <- labels[order.nodes]
  net <- df.to.adj(net)
  net <- net[order.nodes, order.nodes]
  net <- adj.to.df(net)
  net <- net[which(net$weight != 0),]
  
  net$weight <- log10(net$weight) + 7
  
  net$source.label <- labels[net$source]
  net$target.label <- labels[net$target]
  net$source.label <- factor(net$source.label, levels = rev(labels))
  net$target.label <- factor(net$target.label, levels = labels)
  
  pp <- ggplot2::ggplot(net, ggplot2::aes(target.label, source.label, fill=weight))+
    ggplot2::geom_raster(hjust = 0, vjust = 1)+
    viridis::scale_fill_viridis(direction=1, option = 'A')+
    ggplot2::geom_vline(xintercept = bourders, color="#3DED97", size=0.5)+
    ggplot2::geom_hline(yintercept = bourders, color="#3DED97", size=0.5)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size=4),
                   axis.text.y = ggplot2::element_text(size=4))
  
  png(sprintf("%s/%s/%s_linkcommunities.png", path, subfolder, filename), width = 6, height = 5, units = 'in', res = 200)
  print(pp)
  dev.off()
  
  # region.df <- regions[new.idx,]
  # region.df$position <- 1
  # region.df$NODE <- 1:nodes
  # color.regions <- regions[order(regions$REGION),]
  # 
  # gg_leg_y <- ggplot(region.df, aes(y=NODE, x=position, fill=REGION))+
  #   geom_tile()+
  #   scale_x_continuous(expand = c(0, 0))+
  #   scale_y_continuous(trans = 'reverse', expand = c(0, 0))+
  #   scale_fill_manual(values = unique(color.regions$COLOR))+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   theme_void()
  # 
  # leg.gg_y <- get_legend(gg_leg_y)
  # 
  # gg_leg_x <- ggplot(region.df, aes(x=NODE, y=position, fill=REGION))+
  #   geom_tile()+
  #   scale_y_continuous(expand = c(0, 0))+
  #   scale_x_continuous(expand = c(0, 0))+
  #   scale_fill_manual(values = unique(color.regions$COLOR))+
  #   guides(fill=guide_legend(title = 'Area'))+
  #   theme_void()
  # 
  # leg.gg_x <- get_legend(gg_leg_x)
  
  # legend <- plot_grid(leg.pp, leg.gg_y, nrow = 2, rel_heights = c(1,5), axis = 't')
  # legend <- leg.pp
  
  # g <- plot_grid(gg_leg_x + theme(legend.position = 'none',
  #                                 panel.grid = element_blank(),
  #                                 panel.border = element_blank()), 
  #                NULL, NULL,
  #                pp + theme(legend.position = 'none',
  #                           panel.grid = element_blank(),
  #                           panel.border = element_blank()),
  #                gg_leg_y + theme(legend.position = 'none',
  #                                 panel.grid = element_blank(),
  #                                 panel.border = element_blank()), 
  #                legend, 
  #                ncol=3, nrow = 2, rel_widths=c(15,1,2), rel_heights = c(1,15),
  #                align = 'hv' ,
  #                scale=0.98)
}