stan_categories <- function(
  weights, distances, K, serie, inst,
  subfolder="") {
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
        min(),
      weights %>%
        max()
    )
  )
  stan_sample <- rstan::stan(
      file = "STAN/dispersion_problem.stan",
      model_name = "dispersion",
      data = stan_parameters,
      iter = 5000,
      control = list(
        adapt_delta = 0.97,
        max_treedepth = 15
      ),
      verbose = F
    )
  print(
      stan_sample, pars = c(
        "intervals", "lp__"
      )
    )
  sampler_params <- rstan::get_sampler_params(stan_sample, inc_warmup = F)
  summary(do.call(rbind, sampler_params), digits = 2) %>%
    print()
  information <- stan_sample %>%
    rstan::extract()
  intervals <- colMeans(information$intervals)
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
  plot_categories_stan(
    weights, distances, categories, serie, inst,
    subfolder = subfolder
  )
  zeros_categories <- rep(0, length(zeros))
  zeros_categories[!zeros] <- categories
  return(zeros_categories)
}

plot_categories_stan <- function(
  weights, distances, categories, serie, inst,
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
    ggplot2::xlab("Distance") +
    ggplot2::ylab("Weight") +
    ggplot2::theme_bw()
  png("%s/%s/Regression/XGBOOST/%s/%i/stan_categories.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 6, height = 5, res = 200, units = "in")
  print(category_plot)
  dev.off()
}