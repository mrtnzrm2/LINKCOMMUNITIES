measure_similarity_dispersion <- function(
  net, nodes, nt, labels, inst, kfold=3
) {
  source("functions/compute_aik.R")
  source("functions/compute_aki.R")
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/split_kfold.R")
  source("functions/setup_datasets.R")
  exact <- torch::torch_ones(nt, nt, 2)
  # Total error ----
  Error <- c()
  Quad <- c()
  # Start iteration ----
  for (i in 1:100) {
    split_kfold <- split.dataset.kfold(nt, kfold = kfold)
    for (j in 1:kfold) {
      print("** Trial: %i kfold: %i" %>%
        sprintf(i, j)
      )
      ## Get estimate ----
      split <- assemble.split(split_kfold, j)
      datasets <- get.dataset.sim_dist(
        net, nt, nodes, labels, split, inst, zero = T
      )
      sim_df <- datasets$sim_df
      ### Get train sim matrix nodes x ntrn ----
      train_sim <- sim_df$train %>%
        dplyr::as_tibble()
      ntrn <- ncol(train_sim)
      ### Get test sim matrix nodes x (nt - ntrn) ----
      test_sim <- sim_df$test %>%
        dplyr::as_tibble()
      ### Create estimation ----
      estimate <- train_sim %>%
        dplyr::bind_cols(test_sim) %>%
        as.matrix()
      estimate <- estimate[1:nt, 1:nt]
      ## Get exat sim for the k-fold cv ----
      ### Compute exact similarity ----
      net_train <- assemble.dataset.train(
        net %>%
        df.to.adj(),
        split$train, split$test, labels, nt
      )
      net_test <- assemble.dataset.test(
        net %>%
        df.to.adj(),
        split$train, split$test, labels, nt
      )
      net_train <- net_train %>%
        dplyr::bind_cols(net_test) %>%
        as.matrix() %>%
        adj.to.df()
      net_aik <- net_train %>%
        compute.aik(nodes)
      net_aki <- net_train %>%
        compute.aki(nt)
      ## Exact average sim nt x nt ----
      exact[, , 1] <- net_aik[1:nt, 1:nt]
      exact[, , 2] <- net_aki[1:nt, 1:nt]
      exact_mean <- exact %>%
        apply(
          c(1, 2),
          mean
        )
      ## Get quadrants ----
      first_q <- matrix(1, nrow = ntrn, ncol = ntrn) %>%
        dplyr::as_tibble()
      seco_q <- matrix(2, nrow = ntrn, ncol = nt - ntrn) %>%
        dplyr::as_tibble()
      thir_q <- matrix(2, nrow = nt - ntrn, ncol = ntrn) %>%
        dplyr::as_tibble()
      fort_q <- matrix(3, nrow = nt - ntrn, ncol = nt - ntrn) %>%
        dplyr::as_tibble()
      first_q <- first_q %>%
        dplyr::bind_cols(seco_q)
      thir_q <- thir_q %>%
        dplyr::bind_cols(fort_q)
      first_q <- first_q %>%
        dplyr::bind_rows(thir_q)
      first_q <- first_q[1:nt, 1:nt] %>%
        as.matrix() %>%
        adj.to.df() %>%
        dplyr::pull(weight)
      # Compare estimation and exact ----
      error <- exact_mean - estimate
      error[diag(nt) == 1] <- NA
      error <- error %>%
        adj.to.df() %>%
        dplyr::pull(weight)
      keep <- c(
        which(error != 0),
        which(!is.na(error))
      ) %>%
        unique()
      ## Append errors ----
      Error <- c(
        Error,
        error[keep]
      )
      Quad <- c(
        Quad,
        first_q[keep]
      )
    }
  }
  # Plot error ----
  tb <- dplyr::tibble(
      error = Error,
      Case = Quad %>% as.factor()
  )
  tb_mean <- tb %>%
    dplyr::group_by(Case) %>%
    dplyr::summarise(
      mean = mean(error) %>%
        round(3),
      sd = sd(error) %>%
        round(3)
    )
  plot_1 <- tb %>%
    ggplot2::ggplot(
      ggplot2::aes(error, fill = Case)
    ) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = ..density..
      ),
      alpha = 0.8,
      color = "gray"
    ) +
    ggplot2::annotation_custom(
      gridExtra::tableGrob(tb_mean),
      xmin = 0.08,
      ymin = 20
    ) +
    ggplot2::xlab("Error") +
    ggplot2::theme_classic()
  plot_2 <- tb %>%
    ggplot2::ggplot(
      ggplot2::aes(error, fill = Case)
    ) +
    ggplot2::geom_density(
      alpha = 0.6
    ) +
    ggplot2::xlab("Error") +
    ggplot2::theme_classic()
  plot <- cowplot::plot_grid(
    plot_1 + ggplot2::theme(
      legend.position = "none"
    ),
    plot_2,
    rel_widths = c(1, 1.2), labels = "auto"
  )
  png(
    "%s/%s/%s/dispersion_sim.png" %>%
      sprintf(inst$plot, inst$folder, "PreAnalysis"),
    width = 13, height = 5, units = "in", res = 200
  )
  print(plot)
  dev.off()
  png(
    "%s/%s/%s/hist_dispersion_sim.png" %>%
      sprintf(inst$plot, inst$folder, "PreAnalysis"),
    width = 6, height = 5, units = "in", res = 200
  )
  print(plot_1)
  dev.off()
}