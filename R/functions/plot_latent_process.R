plot_latent_process <- function(dataset, serie, inst, subfolder="") {
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>%
   sprintf(inst$plot, inst$folder, subfolder, serie),
   showWarnings = F)
  uncertainty <- dataset$uncertainty
  dataset <- dataset$data
  train_data <- rsample::training(dataset)
  test_data <- rsample::testing(dataset)
  train_data$class <- "train"
  test_data$class <- "test"
  train_data$uncertainty <- uncertainty$train
  test_data$uncertainty <- uncertainty$test
  data <- train_data %>%
    dplyr::bind_rows(test_data)
  # 2D scatter plots ----
  simdist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, sim, color = w)
    ) +
    ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5,
      alpha = 0.5
    ) +
    ggpubr::stat_cor(
      ggplot2::aes(
        label = paste(..r.label.., ..rr.label.., sep = "~`,`~")
      ),
      label.x.npc = 0.5
    ) +
    viridis::scale_color_viridis(name = "log10(FLN) + 7") +
    ggplot2::xlab("Standardized latent distance") +
    ggplot2::ylab("Standardized average similarity") +
    ggplot2::ggtitle("Similarity-Distance") +
    ggplot2::theme_classic()
  wdist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w, color = uncertainty)
    ) +
    ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggpubr::stat_cor(
      ggplot2::aes(
        label = paste(..r.label.., ..rr.label.., sep = "~`,`~")
      ),
      label.x.npc = 0.6
    ) +
    viridis::scale_color_viridis(name = "Uncertainty") +
    ggplot2::xlab("Standardized latent distance") +
    ggplot2::ylab("log10(FLN)+7") +
    ggplot2::ggtitle("Weight-Distance") +
    ggplot2::theme_classic()
  wsim <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(sim, w)
    ) +
     ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggpubr::stat_cor(
      ggplot2::aes(
        label = paste(..r.label.., ..rr.label.., sep = "~`,`~")
      )
    ) +
    ggplot2::xlab("Standardized average similarity") +
    ggplot2::ylab("log10(FLN)+7") +
    ggplot2::ggtitle("Weight-Simlarity") +
    ggplot2::theme_classic()
  p <- cowplot::plot_grid(
    simdist, wsim, wdist, labels = "AUTO",
    ncol = 3
  )
  png("%s/%s/Regression/XGBOOST/%s/%i/scatter_latent.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 17, height = 10, res = 200, units = "in")
  print(p)
  dev.off()
  #*** 3D scatter v
  print("*** Print plotly")
  fig <- data %>%
    plotly::plot_ly(
      x = ~sim, y = ~dist, z = ~w, symbol = ~class,
      symbols = c("diamond", "x"), marker = list(size = 3)
    ) %>%
    plotly::add_markers() %>%
    plotly::layout(
      scene = list(xaxis = list(title = "Standardized average similarity"),
      yaxis = list(title = "Latent distance"),
      zaxis = list(title = "Weight"))
    )
  htmlwidgets::saveWidget(
    fig,
    "%s/%s/Regression/XGBOOST/%s/%i/scatter_latent_3D.html" %>%
      sprintf(inst$plot, inst$folder, subfolder, serie)
  )
  #*** Latent distance vs. weights
  ldist_w <- data %>%
    dplyr::filter(w > 0) %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w, color = class)
    ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.5) +
    ggplot2::geom_smooth(
      method = "lm"
    ) +
    ggpubr::stat_cor() +
    ggplot2::theme_bw()
  png("%s/%s/Regression/XGBOOST/%s/%i/scatter_wldist.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 6, height = 5, res = 200, units = "in")
  print(ldist_w)
  dev.off()
}