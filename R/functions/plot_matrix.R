plot_matrix <- function(M, name1=''){
  source('functions/adj_to_df.R')
  M <- adj.to.df(M)
  M$type <- name1
  p <- ggplot2::ggplot(M, ggplot2::aes(target, source, fill=weight))+
    viridis::scale_fill_viridis(name="prob",option ="A", na.value="white")+
    ggplot2::geom_raster()+
    ggplot2::scale_y_continuous(trans = "reverse")+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  return(p)
}