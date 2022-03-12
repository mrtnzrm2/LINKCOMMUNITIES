plot.aki.aik.density<- function(AKI, AIK, labels, regions, path="", foldername="", subfolder="", suffix="", on=T){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  ## Convert AKI and AIK into data-frames with labels and ordering the labels 
  ## by regions
  source("functions/adj_to_df.R")
  AIK <- adj.to.df(AIK)
  AIK$source.labels <- labels[AIK$source]
  AIK$target.labels <- labels[AIK$target]
  AIK$target.labels <- factor(AIK$target.labels, levels = regions$AREA)
  AIK$source.labels <- factor(AIK$source.labels, levels = rev(regions$AREA))
  
  AKI <- adj.to.df(AKI)
  AKI$source.labels <- labels[AKI$source]
  AKI$target.labels <- labels[AKI$target]
  AKI$target.labels <- factor(AKI$target.labels, levels = regions$AREA)
  AKI$source.labels <- factor(AKI$source.labels, levels = rev(regions$AREA))
  
  AIK$type <- "aik"
  AKI$type <- "aki"
  A <- AIK %>% rbind(AKI)
  
  ph <- ggplot2::ggplot(A, ggplot2::aes(weight, fill=type))+
    ggplot2::geom_density(alpha=0.5, na.rm = T)+
    ggplot2::xlab("similarity")+
    ggplot2::theme_classic()
  
  png(sprintf("%s/%s/%s/density%s.png", path, foldername, subfolder, suffix), width = 6, height = 5, units = 'in', res = 200)
  print(ph)
  dev.off()
}
