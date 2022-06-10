theme_wsim <- function() {
  ggplot2::theme(
    legend.position = "none",
    plot.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#E7DDD7"
    ),
    panel.background = ggplot2::element_rect(
      fill = "#E7DDD7"
    ),
    axis.title = ggplot2::element_text(
      color = "#C00000"
    ),
    axis.text = ggplot2::element_text(
      color = "#4472C4"
    ),
    axis.line = ggplot2::element_line(
      color = "#4472C4"
    ),
    axis.ticks = ggplot2::element_line(
      color = "#4472C4"
    ),
    strip.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#4472C4"
    ),
    strip.text = ggplot2::element_text(
      color = "#C00000"
    ),
    plot.title = ggplot2::element_text(
      color = "#4472C4"
    )
  )
}

theme_b <- function() {
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#E7DDD7"
    ),
    panel.background = ggplot2::element_rect(
      fill = "#E7DDD7"
    ),
    axis.title = ggplot2::element_text(
      color = "#C00000"
    ),
    axis.text = ggplot2::element_text(
      color = "#4472C4"
    ),
    axis.line = ggplot2::element_line(
      color = "#4472C4"
    ),
    axis.ticks = ggplot2::element_line(
      color = "#4472C4"
    ),
    legend.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#E7DDD7"
    ),
    legend.text = ggplot2::element_text(
      color = "#4472C4"
    ),
    legend.title = ggplot2::element_text(
      color = "#C00000"
    ),
    strip.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#4472C4"
    ),
    strip.text = ggplot2::element_text(
      color = "#C00000"
    ),
    plot.title = ggplot2::element_text(
      color = "#4472C4"
    )
  )
}

plot_latent_process_b <- function(dataset, serie, inst, subfolder="") {
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
      color = "#C00000",
      label.x.npc = 0.5
    ) +
    viridis::scale_color_viridis(
      name = "log10(FLN) + 7",
      option = "A"
    ) +
    ggplot2::xlab("Standardized latent distance") +
    ggplot2::ylab("Standardized average similarity") +
    ggplot2::ggtitle("Similarity-Distance") +
    ggplot2::theme_classic() +
    theme_b()
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
      label.x.npc = 0.6,
      color = "#C00000"
    ) +
    viridis::scale_color_viridis(
      name = "Uncertainty",
      option = "A"
    ) +
    ggplot2::xlab("Standardized latent distance") +
    ggplot2::ylab("log10(FLN)+7") +
    ggplot2::ggtitle("Weight-Distance") +
    ggplot2::theme_classic() +
    theme_b()
  wsim <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(sim, w)
    ) +
     ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5,
      alpha = 0.5,
      color = "#4472C4"
    ) +
    ggpubr::stat_cor(
      ggplot2::aes(
        label = paste(..r.label.., ..rr.label.., sep = "~`,`~")
      ),
      color = "#C00000"
    ) +
    ggplot2::xlab("Standardized average similarity") +
    ggplot2::ylab("log10(FLN)+7") +
    ggplot2::ggtitle("Weight-Simlarity") +
    ggplot2::theme_classic() +
    theme_wsim()
  p <- cowplot::plot_grid(
    simdist, wsim, wdist, labels = "AUTO",
    ncol = 3,
    label_colour = "#C00000"
  ) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(
      fill = "#E7DDD7",
      color = "#E7DDD7"
    )
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