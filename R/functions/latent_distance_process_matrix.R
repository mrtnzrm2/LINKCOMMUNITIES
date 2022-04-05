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
      v[v != 0],
      na.rm = T, probs = pracma::linspace(1 / k, 1, k)
    )
  )[1:(k - 1)]
  p[k + 1] <- Inf
  v[v == 0] <- NA
  v <- transform_list_inverse(v, p, length(v)) + 1
  v[is.na(v)] <- 0
  return(v)
}

colSd <- function(x)
  sqrt(rowMeans((t(x) - colMeans(x))^2) * ((dim(x)[1]) / (dim(x)[1] - 1)))

latent_distance_process_matrix <- function(data, st, serie, inst, save=T) {
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/stan_categories.R")
  # Values ----
  categories <- 5
  latent_dimensions <- 3
  nodes_network <- 107
  # Separate data ----
  ## Find categories
  ## Train ----
  train_data <- st$train %>%
    with(
      dplyr::tibble(
        source = source,
        target = target,
        weight = rsample::training(data) %>%
          dplyr::pull(w) %>%
          stan_categories(
            rsample::training(data) %>%
              dplyr::pull(dist),
            categories, serie, inst,
            subfolder = "VariableEvolution")
          )
    ) %>%
    df.to.adj()
  train_dist <- st$train %>%
    with(
      dplyr::tibble(
        source = source,
        target = target,
        weight = rsample::training(data) %>%
            dplyr::pull(dist)
      )
    ) %>%
    df.to.adj()
  col_train <- dim(train_data)[2]
  ## Test ----
  test_data <- st$test %>%
    with(
      dplyr::tibble(
        source = source,
        target = target,
        weight = 0
      )
    ) %>%
    df.to.adj()
  test_dist <- st$test %>%
    with(
      dplyr::tibble(
        source = source,
        target = target,
        weight = rsample::testing(data) %>%
            dplyr::pull(dist)
      )
    ) %>%
    df.to.adj()
  col_test <- dim(test_data)[2]
  ## Create matrice ----
  data_matrix <- matrix(0, nrow = nodes_network, ncol = col_train + col_test)
  data_matrix[, 1:col_train] <- train_data
  data_matrix[, (col_train + 1):(col_train + col_test)] <- test_data
  ## Non-links ----
  no_links <- data_matrix == 0 %>%
    as.integer()
  data_matrix[data_matrix == 0] <- max(data_matrix)
  dist_matrix <- matrix(0, nrow = nodes_network, ncol = col_train + col_test)
  dist_matrix[, 1:col_train] <- train_dist
  dist_matrix[, (col_train + 1):(col_train + col_test)] <- test_dist
  col_data <-  dim(data_matrix)[2]
  if (save) {
    # Format data ----
    stan_set <- list(
      K = categories,
      d = latent_dimensions,
      M = nodes_network,
      N = col_data,
      mat = data_matrix,
      w = no_links,
      D = dist_matrix
    )
    stan_sample <- rstan::stan(
      file = "STAN/zz_org_old.stan",
      model_name = "zz_org_model_old",
      data = stan_set,
      iter = 500,
      control = list(
        adapt_delta = 0.94,
        max_treedepth = 8
      ),
      verbose = F
    )
    print(
      stan_sample, pars = c(
        "sigma", "rho", "rho_s",
        "b", "lp__", "asym[1,2]", "asym[2,1]"
      )
    )
    sampler_params <- rstan::get_sampler_params(stan_sample, inc_warmup = F)
    summary(do.call(rbind, sampler_params), digits = 2) %>%
      print()
    information <- stan_sample %>%
      rstan::extract()
    saveRDS(
      information,
      "../RDS/imputation/original/normal/LatentDistance/4.rds"
    )
    sigma <- information$sigma %>%
      mean()
    asym <- information$asym %>%
      apply(c(2, 3), mean)
    asym_sd <- information$asym %>%
      apply(c(2, 3), sd)
    ldist <- (dist_matrix + asym) / sigma
    train_data <- rsample::training(data)
    train_data$dist <- ldist[, 1:col_train] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    test_data <- rsample::testing(data)
    test_data$dist <- ldist[, (col_train + 1):col_data] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    uncertainty_train <- asym_sd[, 1:col_train] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    uncertainty_test <- asym_sd[, (col_train + 1):col_data] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    ## Transform to data frame and filter zeros in train set ----
    zeros <- train_data$w == 0
    train_data <- train_data %>%
      dplyr::filter(w > 0)
    uncertainty_train <-  uncertainty_train[!zeros]
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
      "../RDS/imputation/original/normal/LatentDistance/4.rds"
    )
    col_data <- 50
    sigma <- information$sigma %>%
      mean()
    asym <- information$asym %>%
      apply(c(2, 3), mean)
    asym_sd <- information$asym %>%
      apply(c(2, 3), sd)
    ldist <- (dist_matrix + asym) / sigma
    train_data <- rsample::training(data)
    train_data$dist <- ldist[, 1:col_train] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    test_data <- rsample::testing(data)
    test_data$dist <- ldist[, (col_train + 1):col_data] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    uncertainty_train <- asym_sd[, 1:col_train] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    uncertainty_test <- asym_sd[, (col_train + 1):col_data] %>%
      adj.to.df() %>%
      dplyr::pull(weight)
    ## Transform to data frame and filter zeros in train set ----
    zeros <- train_data$w == 0
    train_data <- train_data %>%
      dplyr::filter(w > 0)
    uncertainty_train <-  uncertainty_train[!zeros]
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