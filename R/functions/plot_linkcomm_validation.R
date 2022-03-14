plot.linkcommm.validation <- function(train, test, K.coords, serie, inst, subfolder=""){
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>% sprintf(inst$plot, inst$folder, subfolder, serie), showWarnings = F)
  K.coords <- K.coords %>% dplyr::as_tibble()
  colnames(K.coords) <- c("dist", "sim")
  
  p.simdist <- ggplot2::ggplot(train, ggplot2::aes(dist, sim, color=id %>% as.factor()))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::geom_point(data = K.coords, ggplot2::aes(dist, sim), color="black", size=2)+
    ggplot2::ggtitle("Train: similarity-distance")+
    ggplot2::theme_bw()
  
  p.id1 <- train %>% dplyr::select(w, id_1)
  colnames(p.id1) <- c("w", "dist")
  p.id1$id <- "id_1"
  p.id2 <- train %>% dplyr::select(w, id_2)
  colnames(p.id2) <- c("w", "dist")
  p.id2$id <- "id_2"
  # p.id3 <- train %>% dplyr::select(w, id_3)
  # colnames(p.id3) <- c("w", "dist")
  # p.id3$id <- "id_3"
  p <- p.id1 %>% dplyr::bind_rows(p.id2)
  p$cat <- "train"
  
  t.id1 <- test %>% dplyr::select(w, id_1)
  colnames(t.id1) <- c("w", "dist")
  t.id1$id <- "id_1"
  t.id2 <- test %>% dplyr::select(w, id_2)
  colnames(t.id2) <- c("w", "dist")
  t.id2$id <- "id_2"
  # t.id3 <- test %>% dplyr::select(w, id_3)
  # colnames(t.id3) <- c("w", "dist")
  # t.id3$id <- "id_3"
  t <- t.id1 %>% dplyr::bind_rows(t.id2)
  t$cat <- "test"
  
  pt <- p %>% dplyr::bind_rows(t)
  
  p.wid <- ggplot2::ggplot(pt, ggplot2::aes(dist, w, color=cat))+
    ggplot2::facet_wrap(~id)+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::ggtitle("Train: W-centroids distance")+
    ggplot2::theme_bw()
  p.simdist <- cowplot::plot_grid(p.simdist, NULL, ncol = 2, rel_widths = c(3,2))
  p <- cowplot::plot_grid(p.simdist, p.wid, nrow = 2, labels = "auto")
  png("%s/%s/Regression/XGBOOST/%s/%i/wg_test_id.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie), width = 10, height =10 , res = 200, units = "in")
  print(p)
  dev.off()
}