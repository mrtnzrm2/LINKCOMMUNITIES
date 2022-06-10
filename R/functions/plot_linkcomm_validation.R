plot.linkcommm.validation <- function(
  train, test, K.coords, K, serie, inst,
  subfolder="") {
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>%
   sprintf(inst$plot, inst$folder, subfolder, serie),
   showWarnings = F)
  K.coords <- K.coords %>%
    dplyr::as_tibble()
  colnames(K.coords) <- c("dist", "sim")
  train$class <- "train"
  test$class <- "test"
  data <- train %>%
    dplyr::bind_rows(test)
  # Scatter plot of dist and sim with link communities and centroids ----
  simdist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, sim, color = id %>% as.factor())
    ) +
    ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = K.coords,
      ggplot2::aes(dist, sim), color = "black", size = 2
    ) +
    ggplot2::xlab("Standardized distance") +
    ggplot2::ylab("Standardized average similarity") +
    ggplot2::ggtitle("Similarity-Distance") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none"
    )
  wdist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w, color = id %>% as.factor())
    ) +
    ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggplot2::xlab("Standardized distance") +
    ggplot2::ylab("Weight") +
    ggplot2::ggtitle("Weight-Distance") +
    ggplot2::theme_bw()
  wsim <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(sim, w, color = id %>% as.factor())
    ) +
     ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggplot2::xlab("Standardized average similarity") +
    ggplot2::ylab("Weight") +
    ggplot2::ggtitle("Weight-Simlarity") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none"
    )
  p <- cowplot::plot_grid(
    simdist, wsim, wdist,
    ncol = 3, rel_widths = c(1.5, 1.5, 2)
  )
  # # Scatter plot of dist-id and w for different ids ----
  # col.train <- colnames(train)
  # col.train <- col.train[grepl("id_", col.train, fixed = T)]
  # train.ids <- dplyr::tibble()
  # for (ids in col.train) {
  #   p.id <- train %>%
  #     dplyr::select(w, {{ids}})
  #   colnames(p.id) <- c("w", "dist")
  #   p.id$id <- ids
  #   train.ids <- train.ids %>%
  #     dplyr::bind_rows(p.id)
  # }
  # train.ids$cat <- "train"
  # test.ids <- dplyr::tibble()
  # for (ids in col.train) {
  #   p.id <- test %>%
  #     dplyr::select(w, {{ids}})
  #   colnames(p.id) <- c("w", "dist")
  #   p.id$id <- ids
  #   test.ids <- test.ids %>%
  #     dplyr::bind_rows(p.id)
  # }
  # test.ids$cat <- "test"
  # datasets.ids <- train.ids %>%
  #   dplyr::bind_rows(test.ids)
  # datasets.ids$id <- factor(
  #   datasets.ids$id,
  #   levels = paste("id", 1:K, sep = "_")
  # )
  # p.wid <- datasets.ids %>%
  #   ggplot2::ggplot(ggplot2::aes(dist, w, color = cat)) +
  #   ggplot2::facet_wrap(~id) +
  #   ggplot2::geom_point(size = 0.5, alpha = 0.5) +
  #   ggplot2::ggtitle("W-centroids distance") +
  #   ggplot2::theme_bw()
  # p.simdist <- cowplot::plot_grid(
  #     p.simdist, NULL, ncol = 2, rel_widths = c(3, 2)
  # )
  # p <- cowplot::plot_grid(
  #   p.simdist, p.wid,
  #   nrow = 2, labels = "auto", rel_heights = c(2, 4)
  # )
  png("%s/%s/Regression/XGBOOST/%s/%i/centroids_id.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 15, height = 10, res = 200, units = "in")
  print(p)
  dev.off()
  #*** 3D scatter v
  print("*** Print plotly")
  fig <- data %>%
    plotly::plot_ly(
      x = ~sim, y = ~dist, z = ~w,
      color = ~as.factor(id), symbol = ~class,
      symbols = c("diamond", "x"), marker = list(size = 3)
    ) %>%
    plotly::add_markers() %>%
    plotly::layout(
      scene = list(xaxis = list(title = "Standardized average similarity"),
      yaxis = list(title = "Standardized distance"),
      zaxis = list(title = "Weight"))
    )
  htmlwidgets::saveWidget(
    fig,
    "%s/%s/Regression/XGBOOST/%s/%i/centroids_3D.html" %>%
      sprintf(inst$plot, inst$folder, subfolder, serie)
  )
}