plot.aks.matrix <- function(AKS, bourder, colors, labels, regions, path="", foldername="", subfolder="", preffix="", suffix=""){
  
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  ## Convert AKI and AIK into data-frames with labels and ordering the labels 
  ## by regions
  source("functions/adj_to_df.R")
  AKS <- adj.to.df(AKS)
  AKS$source.labels <- labels[AKS$source]
  AKS$target.labels <- labels[AKS$target]
  AKS$target.labels <- factor(AKS$target.labels, levels = regions$AREA)
  AKS$source.labels <- factor(AKS$source.labels, levels = rev(regions$AREA))
  
  p <- ggplot2::ggplot(AKS, ggplot2::aes(target.labels, source.labels, fill=weight))+
    ggplot2::geom_raster(na.rm=T, hjust = 1, vjust = 0)+
    viridis::scale_fill_viridis(option = "A", direction = 1, na.value="white")+
    ggplot2::geom_hline(yintercept = bourder, color="green")+
    ggplot2::geom_vline(xintercept = bourder, color="green")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size = 4, color=colors),
                   axis.text.y = ggplot2::element_text(size=4, color=colors %>% rev()))
  
  png(sprintf("%s/%s/%s/%s_matrix%s.png", path, foldername, subfolder, preffix, suffix), width = 6, height = 5, units = 'in', res = 200)
  print(p)
  dev.off()
}