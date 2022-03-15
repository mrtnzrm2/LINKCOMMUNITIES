plot.linkcommm.validation <- function(
  train, test, K.coords, K, serie, inst,
  subfolder="") {
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>%
   sprintf(inst$plot, inst$folder, subfolder, serie),
   showWarnings = F)
  K.coords <- K.coords %>% dplyr::as_tibble()
  colnames(K.coords) <- c("dist", "sim")
  p.simdist <- ggplot2::ggplot(train,
    ggplot2::aes(dist, sim, color = id %>% as.factor())) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::geom_point(data = K.coords,
    ggplot2::aes(dist, sim), color = "black", size = 2) +
    ggplot2::ggtitle("Train: similarity-distance") +
    ggplot2::theme_bw()
  col.train <- colnames(train)
  col.train <- col.train[grepl("id_", col.train, fixed = T)]
  train.ids <- dplyr::tibble()
  for (ids in col.train) {
    p.id <- train %>% dplyr::select(w, {{ids}})
    colnames(p.id) <- c("w", "dist")
    p.id$id <- ids
    train.ids <- train.ids %>% dplyr::bind_rows(p.id)
  }
  train.ids$cat <- "train"
  test.ids <- dplyr::tibble()
  for (ids in col.train) {
    p.id <- test %>% dplyr::select(w, {{ids}})
    colnames(p.id) <- c("w", "dist")
    p.id$id <- ids
    test.ids <- test.ids %>% dplyr::bind_rows(p.id)
  }
  test.ids$cat <- "test"
  datasets.ids <- train.ids %>% dplyr::bind_rows(test.ids)
  datasets.ids$id <- factor(
    datasets.ids$id,
    levels = paste("id", 1:K, sep = "_")
  )
  p.wid <- datasets.ids %>%
    ggplot2::ggplot(ggplot2::aes(dist, w, color = cat)) +
    ggplot2::facet_wrap(~id) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::ggtitle("Train: W-centroids distance") +
    ggplot2::theme_bw()
  p.simdist <- cowplot::plot_grid(
      p.simdist, NULL, ncol = 2, rel_widths = c(3, 2)
  )
  p <- cowplot::plot_grid(
    p.simdist, p.wid,
    nrow = 2, labels = "auto", rel_heights = c(2, 4)
  )
  png("%s/%s/Regression/XGBOOST/%s/%i/wg_test_id.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 10, height = 10, res = 200, units = "in")
  print(p)
  dev.off()
}