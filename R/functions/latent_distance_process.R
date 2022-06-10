transform_list_inverse <- function(A, p, N) {
  K <- length(p)
  for (i in 1:N) {
    if (!is.na(A[i])) {
        for (k in 1:(K - 1)) {
            if (A[i] >= p[k] && A[i] < p[k + 1]) {
            A[i] <- K - 1 - k
            break
            }
        }
    }
  }
  return(A)
}

categorize <- function(v, k) {
  p <- rep(0, k + 1)
  p[2:k] <- unname(
    quantile(
      v,
      na.rm = T, probs = pracma::linspace(1 / k, 1, k)
    )
  )[1:(k - 1)]
  p[k + 1] <- Inf
  v <- transform_list_inverse(v, p, length(v)) + 1
  return(v)
}

colSd <- function(x)
  sqrt(rowMeans((t(x) - colMeans(x))^2) * ((dim(x)[1]) / (dim(x)[1] - 1)))

latent_distance_process <- function(data, st, serie, inst, save=T) {
  source("functions/stan_categories.R")
  # Values ----
  number_categories <- 5
  number_latent_dimensions <- 1
  number_nodes_network <- 107
  # Separate data ----
  ## Train ----
  train_data <- rsample::training(data)
  ntrn <- train_data %>%
    nrow()
  train_cat <- train_data %>%
    dplyr::pull(w)
  train_dist <- train_data %>%
    dplyr::pull(dist)
  train_cat <- train_cat %>%
    stan_categories(
      train_dist, number_categories, serie, inst,
      subfolder = "VariableEvolution")
  train_sim <- train_data %>%
    dplyr::pull(sim)
  train_source <- st$train %>%
    dplyr::pull(source)
  train_target <- st$train %>%
    dplyr::pull(target)
  ## Test ----
  test_data <- rsample::testing(data)
  ntst <- test_data %>%
    nrow()
  test_dist <- test_data %>%
    dplyr::pull(dist)
  test_sim <- test_data %>%
    dplyr::pull(sim)
  test_source <- st$test %>%
    dplyr::pull(source)
  test_target <- st$test %>%
    dplyr::pull(target)
  if (save) {
    # Format data ----
    stan_set <- list(
      K = number_categories,
      d = number_latent_dimensions,
      M = number_nodes_network,
      category_train = train_cat,
      dist_train = train_dist,
      sim_train = train_sim,
      source_train = train_source,
      target_train = train_target,
      dist_test = test_dist,
      sim_test = test_sim,
      source_test = test_source,
      target_test = test_source,
      ntrn = ntrn,
      ntst = ntst
    )
    stan_sample <- rstan::stan(
      file = "STAN/zz_org_2.stan",
      model_name = "zz_org_model_2",
      data = stan_set,
      iter = 500,
      control = list(
        adapt_delta = 0.96,
        max_treedepth = 8
      ),
      verbose = T
    )
    print(
      stan_sample, pars = c(
        "sigma", "rho", "rho_s",
        "b", "lp__", "train_ldist[1]", "test_ldist[2]"
      )
    )
    sampler_params <- rstan::get_sampler_params(stan_sample, inc_warmup = F)
    summary(do.call(rbind, sampler_params), digits = 2) %>%
      print()
    information <- stan_sample %>%
      rstan::extract()
    # saveRDS(
    #   information,
    #   "../RDS/imputation/original/normal/LatentDistance/3.rds"
    # )
    train_data$dist <- information$train_ldist %>%
      colMeans()
    test_data$dist <- information$test_ldist %>%
      colMeans()
    uncertainty_train <- information$train_ldist %>%
      colSd()
    uncertainty_test <- information$test_ldist %>%
      colSd()
    data <- structure(
      list(
        data = train_data %>%
          rbind(test_data) %>%
          dplyr::as_tibble(),
        in_id = nrow(train_data) %>%
          seq_len(),
        out_id = (nrow(train_data) + 1):(nrow(train_data) + nrow(test_data))
      ),
      class = "rsplit"
    )
  } else{
    information <- readRDS(
      "../RDS/imputation/original/normal/LatentDistance/3.rds"
    )
    train_data$dist <- information$train_ldist %>%
      colMeans()
    test_data$dist <- information$test_ldist %>%
      colMeans()
    uncertainty_train <- information$train_ldist %>%
      colSd()
    uncertainty_test <- information$test_ldist %>%
      colSd()
    data <- structure(
      list(
        data = train_data %>%
          rbind(test_data) %>%
          dplyr::as_tibble(),
        in_id = nrow(train_data) %>%
          seq_len(),
        out_id = (nrow(train_data) + 1):(nrow(train_data) + nrow(test_data))
      ),
      class = "rsplit"
    )
  }
  return(
    list(
      data = data,
      uncertainty = list(
        train = uncertainty_train,
        test = uncertainty_test
      )
    )
  )
}