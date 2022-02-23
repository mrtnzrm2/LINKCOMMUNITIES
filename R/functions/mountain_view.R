mountain.view <- function(ham.merde, labels, net.coord, regions, path='', nodes=40,
                          animal='monkey', foldername='', subfolder='',
                          plt=F, filename=''){
  
  source('functions/in_the_heights.R')
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  
  net.coord <- as.matrix(net.coord)
  
  x.coord <- net.coord[,1]
  y.coord <- net.coord[,2]
  z.coord <- in.the.heights(ham.merde, nodes)[1:nodes]
  
  s <- akima::interp(x.coord, y.coord, z.coord, nx = 300, ny = 300)
  
  regions <- regions[which(regions$AREA %in% labels[1:nodes]),]
  regions <- regions[match(labels[1:nodes], regions$AREA),]
  
  data.3d <- data.frame(x=x.coord,
                        y=y.coord,
                        z=z.coord[1:nodes],
                        REGION=as.factor(regions$REGION))
  
  surface.3d <- data.frame(source=rep(s$x, each=300), target=rep(s$y, 300), weight=pracma::Reshape(t(s$z), 300*300, 1))
  surface.3d <- surface.3d[!is.na(surface.3d$weight),]
  
  pp <-ggplot2::ggplot(surface.3d) +
    ggplot2::geom_tile(ggplot2::aes(x = source, y = target, fill = weight)) +
    ggplot2::geom_contour(ggplot2::aes(x = source, y = target, z = weight), color = "gray", 
                 binwidth=10)+
    viridis::scale_fill_viridis(option = 'A')+
    ggplot2::geom_point(data=data.3d, ggplot2::aes(x=x, y=y, color=REGION), size=10)+
    ggplot2::annotate('text', data.3d$x,
             data.3d$y, label=labels[1:nodes],
             color='white', fontface=2, size=4)+
    ggplot2::coord_fixed()+
    ggplot2::theme_void()
  
  
  if (plt){
    png(sprintf("%s/%s/%s/%s.png", path, foldername, subfolder, filename), width = 13, height = 10, units = 'in', res = 200)
    print(pp)
    dev.off()  
  } else {
    print(pp)
  }
  
}