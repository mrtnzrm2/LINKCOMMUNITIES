plot.parameters <- function(Dnetwork, best.height, view='Dview', filename='', plt = T, animal='monkey',
                            foldername='', subfolder='', path=''){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  dot.plot <- Dnetwork[,c('height',view)]
  colnames(dot.plot) <- c('h', 'w')
  p2 <- ggplot2::ggplot(dot.plot, ggplot2::aes(h,w))+
    ggplot2::geom_point(color='red', size=0.1)+
    ggplot2::geom_line()+
    ggplot2::annotate('text', x=best.height, y=0, label=format(best.height), cex=5)+
    ggplot2::geom_vline(xintercept = best.height, color='blue')+
    ggplot2::ggtitle(filename)+
    ggplot2::ylab(view)+
    ggplot2::theme_classic()
  
  if(plt){
    png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 8, height = 6, units = 'in', res = 200)
    print(p2)
    dev.off()
  } else{
    return(p2)
  }
}