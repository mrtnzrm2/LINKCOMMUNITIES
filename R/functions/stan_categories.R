stan_categories <- function(
  weights, distances, K, serie, inst,
  subfolder="", save=T) {
  zeros <- weights == 0
  weights <- weights[!zeros]
  distances <- distances[!zeros]
  order_weights <- sort(weights)
  stan_parameters <- list(
    w = weights[order_weights],
    dist = distances[order_weights],
    n = weights %>%
      length(),
    K = K,
    range = c(
      weights %>%
        min() - 1e-5,
      weights %>%
        max() + 1e-5
    )
  )
  # CmdStanR ----
  ## Fitting ----
  stan_model <-  cmdstanr::cmdstan_model("STAN/dispersion_problem.stan")
  model_MCMC <- stan_model$sample(
      data = stan_parameters,
      parallel_chains = 4,
      adapt_delta = 0.99998,
      max_treedepth = 10,
      iter_warmup = 1000,
      iter_sampling = 1000,
  )
  ## Diagnose ----
  print(
      model_MCMC$summary(
        variables = c(
          "intervals", "lp__"
        )
      )
  )
  print(
      model_MCMC$diagnostic_summary()
  )
  print(
    model_MCMC$cmdstan_diagnose()
  )
  ## Mean and Sd ----
  information <- model_MCMC$summary(
      variables = "intervals",
      "mean"
  )
  intervals <- information %>%
    dplyr::pull(mean)
  Intervals <- c(
    weights %>%
      min(),
    intervals,
    weights %>%
      max()
  )
  categories <- weights %>%
    lapply(
      function(x) {
        for (i in 1:K) {
          if (
            x >= Intervals[i] &&
            x <= Intervals[i + 1]
          )
            return(K - i + 1)
        }
      }
    ) %>%
    unlist()
  if (save)
    plot_categories_stan(
      weights, distances, intervals, categories, serie, inst,
      subfolder = subfolder
    )
  no_zeros_categories <- rep(0, length(zeros))
  no_zeros_categories[!zeros] <- categories
  return(no_zeros_categories)
}

plot_categories_stan <- function(
  weights, distances, intervals, categories, serie, inst,
  subfolder="") {
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>%
   sprintf(inst$plot, inst$folder, subfolder, serie),
   showWarnings = F
  )
  category_plot <- dplyr::tibble(
    w = weights,
    dist = distances,
    category = categories %>%
      as.factor()
  ) %>%
    ggplot2::ggplot(
      ggplot2::aes(dist, w, color = category)
    ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = intervals) +
    ggplot2::xlab("tracto-dist/max(tracto-dist)") +
    ggplot2::ylab("log10(FLN)+7") +
    ggplot2::theme_classic()
  png("%s/%s/Regression/XGBOOST/%s/%i/stan_categories.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 6, height = 5, res = 200, units = "in")
  print(category_plot)
  dev.off()
}