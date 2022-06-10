colSd <- function(x)
  sqrt(rowMeans((t(x) - colMeans(x))^2) * ((dim(x)[1]) / (dim(x)[1] - 1)))

with_naive_categories <- function(data, st, categories, serie, inst, save=T) {
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  # weight <- rsample::training(data)
  # print(st$train %>% head())
  # print(weight %>% head())
  weight <-  rsample::training(data) %>%
    dplyr::pull(w) %>%
    naive_categories(
      rsample::training(data) %>%
        dplyr::pull(dist),
      categories, serie, inst,
      subfolder = "VariableEvolution",
      save = save
    )
  return(
    st$train %>%
      with(
        dplyr::tibble(
          source = source,
          target = target,
          weight = weight
          )
      ) %>%
      df.to.adj()
  )
}

with_stan_categories <- function(data, st, serie, categories, inst, save=T) {
  source("functions/df_to_adj.R")
  return(
    st$train %>%
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
              subfolder = "VariableEvolution",
              save = save
            )
          )
      ) %>%
      df.to.adj()
  )
}

latent_distance_process_matrix <- function(data, st, serie, inst, save=T) {
  source("functions/df_to_adj.R")
  source("functions/adj_to_df.R")
  source("functions/stan_categories.R")
  source("functions/naive_categories.R")
  # Values ----
  categories <- 5
  latent_dimensions <- 3
  nodes_network <- 107
  # Separate data ----
  ## Find categories
  ## Train ----
  # train_data <- with_stan_categories(
  #   data, st, serie, categories, inst, save = T
  # )
  train_data <- with_naive_categories(
    data, st, categories, serie, inst, save = T
  )
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
  ## Create matrices ----
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
  ### sigma ----
  sigma <- 2
  if (save) {
    while (sigma > 1) {
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
      # CmdstanR ----
      stan_model <- cmdstanr::cmdstan_model("STAN/zz_org_old.stan")
      model_MCMC <- stan_model$sample(
        data = stan_set,
        parallel_chains = 4,
        adapt_delta = 0.95,
        max_treedepth = 9,
        iter_warmup = 250,
        iter_sampling = 250,
      )
      ### CmdstanR summary ----
      print(
        model_MCMC$summary(
          variables = c(
            "sigma", "rho", "rho_s",
            "b", "lp__", "asym[1,2]", "asym[2,1]"
          )
        )
      )
      print(
        model_MCMC$diagnostic_summary()
      )
      print(
        model_MCMC$cmdstan_diagnose()
      )
      ### Mean and Sd -----
      information <- model_MCMC$summary(
        variables = c("sigma", "asym"),
        "mean", "sd"
      )
      saveRDS(
        information,
        "../RDS/imputation/original/normal/LatentDistance/4.rds"
      )
      sigma <- information %>%
        dplyr::filter(variable == "sigma") %>%
        dplyr::pull(mean)
    }
    asym <- information %>%
      dplyr::filter(grepl("asym", variable)) %>%
      dplyr::pull(mean) %>%
      pracma::Reshape(nodes_network, col_data)
    asym_sd <- information %>%
      dplyr::filter(grepl("asym", variable)) %>%
      dplyr::pull(sd) %>%
      pracma::Reshape(nodes_network, col_data)
    source("functions/standardized.R")
    ldist <- (dist_matrix + asym) / sigma
    ldist <- ldist %>%
      standardized.block.diag()
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
    # Manually set up: col_data!!
    col_data <- 50
    sigma <- information$mean[information$variable == "sigma"]
    asym <- information %>%
      dplyr::filter(grepl("asym", variable, fixed = T)) %>%
      dplyr::pull(mean) %>%
      pracma::Reshape(nodes_network, col_data)
    asym_sd <- information %>%
      dplyr::filter(grepl("asym", variable, fixed = T)) %>%
      dplyr::pull(sd) %>%
      pracma::Reshape(nodes_network, col_data)
    source("functions/standardized.R")
    ldist <- (dist_matrix + asym) / sigma
    ldist <-  ldist %>%
      standardized.block.diag()
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