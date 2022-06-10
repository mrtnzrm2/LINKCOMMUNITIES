plot.xgb.residuals <- function(serie, inst, subfolder="", suffix="") {
  source("functions/rmae.R")
  dir.create(
        "%s/%s/Regression/XGBOOST/%s/%i" %>%
            sprintf(inst$plot, inst$folder, subfolder, serie),
        showWarnings = F
    )
  # Load data ----
  data <- readRDS(
      "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
        sprintf(inst$common, serie)
  )
  # Marginalize over ids
  data_1 <- data %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      w = unique(w)[1],
      stdist = unique(dist)[1],
      .rmae = mean(abs(w - pred) / w),
      .pred = mean(pred),
      .rmse = Metrics::rmse(w, pred),
      .sim  = mean(sim)
    ) %>%
    dplyr::mutate(
      .resid = (w - .pred),
      .stdresid = (w - .pred) / sd((w - .pred))
    )
  source("functions/plot_diagnosis.R")
  # Plot diagnosis ----
  p <- plot.diagnosis(data_1)
  purb <- ggpubr::ggarrange(
    p$res.fit, p$norm.qq, p$scale.location, p$den.res,
    nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")
  )
  png(
    "%s/%s/Regression/XGBOOST/%s/%i/residuals_marg_id%s.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie, suffix),
    width = 9, height = 8, res = 200, units = "in"
  )
  print(purb)
  dev.off()
  # Plot RMAE ----
 rmae.table <- dplyr::tibble(
    w = c("Weak", "Medium", "Strong"),
    rmae = c(
              data_1 %>%
                dplyr::filter(w <  2.5) %>%
                dplyr::pull(.rmae) %>%
                mean() %>%
                round(2),
              data_1 %>%
                dplyr::filter(w >= 2.5, w <= 4.5) %>%
                dplyr::pull(.rmae) %>%
                mean() %>%
                round(2),
              data_1 %>%
                dplyr::filter(w > 4.5) %>%
                dplyr::pull(.rmae) %>%
                mean() %>%
                round(2)
    )
  )
  p.wrmae <- ggplot2::ggplot(data_1, ggplot2::aes(w, .rmae)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("w") +
    ggplot2::ylab("Relative prediction error") +
    ggplot2::ylim(-0.1, 3) +
    ggplot2::stat_smooth(method = "loess", color = "orange", fill = "orange") +
    ggplot2::annotation_custom(
      gridExtra::tableGrob(rmae.table),
      xmin = 5, xmax = 5.5, ymin = 2, ymax = 2.25
    ) +
    ggplot2::theme_bw()
  p.distrmae <- ggplot2::ggplot(data_1, ggplot2::aes(stdist, .rmae)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("Standardize distance") +
    ggplot2::ylab("Relative prediction error") +
    ggplot2::stat_smooth(
      method = "loess", color = "orange", fill = "orange", formula = y ~ x
    ) +
    ggplot2::theme_bw()
  p.simrmae <- ggplot2::ggplot(data_1, ggplot2::aes(.sim, .rmae)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("Standardize similarity") +
    ggplot2::ylab("Relative prediction error") +
    ggplot2::stat_smooth(
      method = "loess", color = "orange", fill = "orange", formula = y ~ x
    ) +
    ggplot2::theme_bw()
  p.wrmse <- ggplot2::ggplot(data_1, ggplot2::aes(w, .rmse)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("w" ) +
    ggplot2::ylab("Root mean square error") +
    ggplot2::stat_smooth(
      method = "loess", color = "orange", fill = "orange", formula = y ~ x
    ) +
    ggplot2::theme_bw()
  p.distrmse <- ggplot2::ggplot(data_1, ggplot2::aes(stdist, .rmse)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("Standardize distance") +
    ggplot2::ylab("Root mean square error") +
    ggplot2::stat_smooth(
      method = "loess", color = "orange", fill = "orange", formula = y ~ x
    )+
    ggplot2::theme_bw()
  p.simrmse <- ggplot2::ggplot(data_1, ggplot2::aes(.sim, .rmse)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::xlab("Standardize similarity") +
    ggplot2::ylab("Root mean square error") +
    ggplot2::stat_smooth(
      method = "loess", color = "orange", fill = "orange", formula = y ~ x
    ) +
    ggplot2::theme_bw()
  purb <- cowplot::plot_grid(
    p.wrmae, p.distrmae, p.simrmae, p.wrmse, p.distrmse, p.simrmse,
    nrow = 2, ncol = 3, labels = "auto"
  )
  png(
    "%s/%s/Regression/XGBOOST/%s/%i/rmae_marg_id%s.png"
      %>%
      sprintf(inst$plot, inst$folder, subfolder, serie, suffix),
      width = 15, height = 9, res = 200, units = "in"
  )
  print(purb)
  dev.off()
  # Plot train-test histogram ----
  #*** Get mean score from k-fold v
  histo_train_test <- data %>%
    dplyr::group_by(trial) %>%
    dplyr::summarise(
      train_rmae = mean(train_rmae, na.rm = T),
      test_rmae = mean(test_rmae, na.rm = T)
    )
  #*** Get distance v
  histo_train_test <- histo_train_test %>%
    dplyr::mutate(
      dist = (train_rmae - test_rmae) / test_rmae
    )
  #*** Extract and create tibble witht two categories v
  histo_1 <- dplyr::bind_rows(
    dplyr::tibble(
      weight = histo_train_test %>%
        dplyr::pull(train_rmae),
      type = "train"
    ),
    dplyr::tibble(
      weight = histo_train_test %>%
        dplyr::pull(test_rmae),
      type = "test"
    )
  )
  #*** Get mean scores v
  mean_histo_1 <- histo_1 %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(
      weight = mean(weight, na.rm = T)
    )
  #*** T-test v
  t_test <- t.test(
    histo_train_test %>%
      dplyr::pull(train_rmae),
    histo_train_test %>%
      dplyr::pull(test_rmae),
      var.equal = F
  )
  t_test <- dplyr::tibble(
    weight = t_test$p.value
  )
  #*** Plot histogram v
  p_hist_1 <- histo_1 %>%
    ggplot2::ggplot(
      ggplot2::aes(weight, fill = type)
    ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ..density..),
      position = "identity",
      alpha = 0.5,
      color = "gray",
      bins = 20
    ) +
    ggplot2::geom_vline(
      xintercept = mean_histo_1$weight,
      color = "blue",
      alpha =  0.5,
      linetype = "dashed"
    ) +
    ggplot2::geom_density(
      alpha = 0.3
    ) +
    ggplot2::geom_text(
      data = mean_histo_1,
      ggplot2::aes(
        x = mean_histo_1$weight,
        y = 350,
      ),
      label = mean_histo_1$weight %>%
        round(3),
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = 0.245,
      y = 350,
      color = "red",
      label = "p < %s" %>%
        sprintf(
          t_test$weight %>%
            formatC(
              format = "e",
              digits = 2
          )
        )
    ) +
    ggplot2::xlab("RMAE") +
    ggplot2::theme_bw()
  #*** Plot scatter v
  p_scatter_1 <- data_1 %>%
    ggplot2::ggplot(
      ggplot2::aes(.pred, w)
    ) +
    ggplot2::geom_point(
      size = 0.5,
      alpha = 0.8
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      color = "red"
    ) +
    ggpubr::stat_cor(
      color = "red"
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "blue",
      alpha = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::xlab("prediction") +
    ggplot2::ylab("w") +
    ggplot2::theme_bw()
  #*** Distance histogram v
  p_hist_2 <- histo_train_test %>%
    ggplot2::ggplot(
      ggplot2::aes(dist)
    ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ..density..),
      bins = 20,
      color = "black",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = histo_train_test$dist %>%
        mean(),
      color = "blue",
      alpha = 0.8,
      linetype = "dashed"
    ) +
    ggplot2::xlab("(train-test)/test") +
    ggplot2::theme_bw()
  p_1 <- cowplot::plot_grid(
    p_hist_1, p_hist_2,
    ncol = 2, labels = "auto", rel_widths = c(3, 2)
  )
  p_2 <- cowplot::plot_grid(
    NA, p_scatter_1, NA,
    ncol = 3, labels = c("", "c", ""), rel_widths = c(1, 3, 1)
  )
  p_3 <- cowplot::plot_grid(
    p_1, p_2,
    nrow = 2, rel_heights = c(3, 2)
  )
  png(
    "%s/%s/Regression/XGBOOST/%s/%i/scores_%s.png"
      %>%
      sprintf(inst$plot, inst$folder, subfolder, serie, suffix),
      width = 9, height = 9, res = 200, units = "in"
  )
  print(p_3)
  dev.off()
}