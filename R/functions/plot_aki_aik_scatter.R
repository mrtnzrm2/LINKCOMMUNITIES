plot.aki.aik.scatter <- function(net, AKI, AIK, nt, labels, path="", foldername="", subfolder="", suffix=""){
  if (!pracma::strcmp(subfolder, '')){
    dir.create(sprintf('%s/%s/%s', path, foldername, subfolder), showWarnings = FALSE)
  }
  ## Convert AKI and AIK into data-frames with labels and ordering the labels 
  ## by regions
  source("functions/adj_to_df.R")
  AIK <- adj.to.df(AIK)
  AIK$weight[AIK$target == AIK$source] <- NA
  AKI <- adj.to.df(AKI)
  AKI$weight[AKI$target == AKI$source] <- NA
  
  ## Get distances
  # source("functions/get_tracto2016.R")
  # distances <- get.tracto2016(labels)
  source("functions/get_lv_tracto2016.R")
  distances <- get.lv.tracto2016(labels)
  distances <- distances %>% adj.to.df()
  distances$weight[distances$target == distances$source] <- NA
  
  ## Format net
  source("functions/df_to_adj.R")
  net <- with(net, data.frame(source=source, target=target, weight=weight))
  net <- net %>% df.to.adj()
  net[net == 0] <- NA
  net <- net %>% adj.to.df()
  
  print("AIK-AKI")
  AKS <- data.frame(aik=AIK$weight[AIK$target <= nt & AIK$source <= nt & AIK$source < AIK$target], aki=AKI$weight[AKI$target <= nt & AKI$source < AKI$target])
  print("AIK-DIS")
  AIKDIS <- data.frame(aik=AIK$weight[AIK$target <= nt & AIK$source <= nt & AIK$source < AIK$target], distance=distances$weight[distances$target <= nt & distances$source <= nt & distances$source < distances$target])
  print("AKI-DIS")
  AKIDIS <- data.frame(aki=AKI$weight[AIK$target <= nt & AIK$source <= nt & AIK$source < AIK$target], distance=distances$weight[distances$target <= nt & distances$source <= nt & distances$source < distances$target])
  print("AIK-W")
  AIKW <- data.frame(aik=AIK$weight[AIK$target <= nt], w=net$weight[net$target <= nt])
  print("AKI-W")
  AKIW <- data.frame(aki=AKI$weight[AKI$target <= nt], w=net$weight[net$target <= nt & net$source <= nt])
  print("W-DIS")
  WDIS <- data.frame(distance=distances$weight[distances$target <= nt], w=net$weight[net$target <= nt])
  
  pks <- ggplot2::ggplot(AKS, ggplot2::aes(aik, aki))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red")+
    ggplot2::ggtitle("AKI-AIK")
  paikdis <- ggplot2::ggplot(AIKDIS, ggplot2::aes(aik, distance))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red")+
    ggplot2::ggtitle("AIK-DIS")
  pakidis <- ggplot2::ggplot(AKIDIS, ggplot2::aes(aki, distance))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red")+
    ggplot2::ggtitle("AKI-DIS")
  pakiw <- ggplot2::ggplot(AKIW, ggplot2::aes(aki, w))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red")+
    ggplot2::ggtitle("AKI-W")
  paikw <- ggplot2::ggplot(AIKW, ggplot2::aes(aik, w))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red", label.y = 0.5)+
    ggplot2::ggtitle("AIK-W")
  pwdis <- ggplot2::ggplot(WDIS, ggplot2::aes(distance, w))+
    ggplot2::geom_point(na.rm = T, size=0.5)+
    ggplot2::geom_smooth(method = "lm")+
    ggpubr::stat_cor(color="red", label.x =7.5, label.y = 0.5)+
    ggplot2::ggtitle("W-DIS")
  
  purb <- ggpubr::ggarrange(pks, paikdis, pakidis, paikw, pakiw, pwdis, nrow = 2, ncol = 3, labels = c("A", "B", "C", "D", "E", "F"))
  
  png(sprintf("%s/%s/%s/scatter%s.png", path, foldername, subfolder, suffix), width = 12, height = 6, units = 'in', res = 200)
  print(purb)
  dev.off()
}