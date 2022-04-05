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
    ggpubr::stat_cor() +
    viridis::scale_color_viridis() +
    ggplot2::xlab("Distance") +
    ggplot2::ylab("Standardized average similarity") +
    ggplot2::ggtitle("Similarity-Distance") +
    ggplot2::theme_bw()
  wdist <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w, color = uncertainty)
    ) +
    ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggpubr::stat_cor() +
    viridis::scale_color_viridis(option = "A") +
    ggplot2::xlab("Latent distance") +
    ggplot2::ylab("Weight") +
    ggplot2::ggtitle("Weight-Distance") +
    ggplot2::theme_bw()
  wsim <- data %>%
    ggplot2::ggplot(
      ggplot2::aes(sim, w, color = dist)
    ) +
     ggplot2::facet_wrap(~class, nrow = 2) +
    ggplot2::geom_point(
      size = 0.5, alpha = 0.5
    ) +
    ggpubr::stat_cor() +
    viridis::scale_color_viridis(option = "A") +
    ggplot2::xlab("Standardized average similarity") +
    ggplot2::ylab("Weight") +
    ggplot2::ggtitle("Weight-Simlarity") +
    ggplot2::theme_bw()
  p <- cowplot::plot_grid(
    simdist, wsim, wdist,
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

just_distance <- function(dataset) {
  train <- rsample::training(dataset) %>%
    dplyr::select(-sim)
  test <- rsample::testing(dataset) %>%
    dplyr::select(-sim)
  data <- structure(
    list(
      data = train %>%
        rbind(test) %>%
        dplyr::as_tibble(),
      in_id = nrow(train) %>%
        seq_len(),
      out_id = (nrow(train) + 1):(nrow(train) + nrow(test))
    ),
    class = "rsplit"
  )
  return(data)
}

run_ogr_regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_ogr.R")
    source("functions/latent_distance_process.R")
    source("functions/latent_distance_process_matrix.R")
    print("* Starting xgb regression")
    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    # Start experiment ----
    for (i in 1:100) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        print("** Trial: %i kfold: %i" %>%
          sprintf(i, j)
        )
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset.sim_dist(
          net, nt, nodes, labels, split, inst, zero = T
        )
        ids <- datasets$ids
        st <- datasets$st
        datasets <- datasets$data
        datasets <- latent_distance_process_matrix(
          datasets, st, serie, inst,
          save = T
        )
        plot_latent_process(
          datasets, serie, inst, subfolder = "VariableEvolution"
        )
        datasets <- datasets$data
        ### Run XGBOOST ----
        xgboost_model <- parsnip::boost_tree(
            mode = "regression",
            trees = 750, # 750 1000
            min_n = 18, # 18 20
            tree_depth =  2, # 2 15
            learn_rate = 0.0094047, # 0.0094047 0.0142631
            loss_reduction = 0.0458003 # 0.0458003 0.1355005
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        # ### Save ----
        save_work_ogr(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}

set_subserie <- function(subserie, n=100, m=10) {
  long_serie <- 1:n
  return(
    long_serie[(((subserie - 1) * m) + 1):(subserie * m)]
  )
}

run_ogr_regression_crc <- function(
  net, nt, nodes, labels, serie, subserie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_ogr.R")
    source("functions/latent_distance_process.R")
    source("functions/latent_distance_process_matrix.R")
    print("* Starting xgb regression")
    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    # Start experiment ----
    for (i in set_subserie(subserie)) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        print("** Trial: %i kfold: %i" %>%
          sprintf(i, j)
        )
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset.sim_dist(
          net, nt, nodes, labels, split, inst, zero = T
        )
        ids <- datasets$ids
        st <- datasets$st
        datasets <- datasets$data
        datasets <- latent_distance_process_matrix(
          datasets, st, serie, inst,
          save = T
        )
        plot_latent_process(
          datasets, serie, inst, subfolder = "VariableEvolution"
        )
        datasets <- datasets$data
        ### Run XGBOOST ----
        xgboost_model <- parsnip::boost_tree(
            mode = "regression",
            trees = 750, # 750 1000
            min_n = 18, # 18 20
            tree_depth =  2, # 2 15
            learn_rate = 0.0094047, # 0.0094047 0.0142631
            loss_reduction = 0.0458003 # 0.0458003 0.1355005
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        # ### Save ----
        save_work_ogr(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}