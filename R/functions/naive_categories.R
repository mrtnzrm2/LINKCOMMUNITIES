transform_list_inverse <- function(A, p, N) {
  categories <- rep(0, N)
  K <- length(p)
  for (i in 1:N) {
    if (!is.na(A[i])) {
        for (k in 1:(K - 1)) {
            if (A[i] >= p[k] && A[i] < p[k + 1]) {
            categories[i] <- K - 1 - k
            break
            }
        }
    }
  }
  return(categories)
}

naive_categories <- function(
  weight, dist, k, serie, inst,
  subfolder="", save=T
) {
  intervals <- rep(0, k + 1)
  intervals[2:k] <- unname(
    quantile(
      weight[weight != 0],
      na.rm = T, probs = pracma::linspace(1 / k, 1, k)
    )
  )[1:(k - 1)]
  intervals[k + 1] <- Inf
  weight[weight == 0] <- NA
  categories <- transform_list_inverse(weight, intervals, length(weight)) + 1
  if (save) {
    plot_categories_naive(
      weight, dist, intervals[2:k], categories,
      serie, inst, subfolder = subfolder
    )
  }
  categories[is.na(weight)] <- 0
  return(categories)
}

plot_categories_naive <- function(
  weight, distances, intervals, categories, serie, inst,
  subfolder="") {
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>%
   sprintf(inst$plot, inst$folder, subfolder, serie),
   showWarnings = F
  )
  category_plot <- dplyr::tibble(
    w = weight,
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
  png("%s/%s/Regression/XGBOOST/%s/%i/naive_categories.png" %>%
    sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 6, height = 5, res = 200, units = "in")
  print(category_plot)
  dev.off()
}