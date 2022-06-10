wlm_weight_pred <- function(
  train_w, train_sim, train_dist, test_sim, test_dist
  ) {
  raw_model <- lm(
    w ~ sim + dist,
    data = dplyr::tibble(
      w = train_w,
      sim = train_sim,
      dist = train_dist
    )
  )
  wt <- 1 / lm(
    abs(raw_model$residuals) ~ raw_model$fitted.values
  )$fitted.values^2
  wls_model <- lm(
    w ~ sim + dist,
    data =  dplyr::tibble(
      w = train_w,
      sim = train_sim,
      dist = train_dist
    ),
    weights = wt
  )
  return(
    predict(
      wls_model,
      dplyr::tibble(
        sim = test_sim,
        dist = test_dist
      )
    )
  )
}

link.classification <- function(dataset, zeros, serie, inst, subfolder="") {
  source("functions/df_to_adj.R")
  source("functions/standardized.R")
  # Ploting zeros information ----
  make_zeros(
    rsample::training(dataset), rsample::testing(dataset), zeros$K_coords,
    serie, inst,
    subfolder = subfolder, on = F
  )
  # Take out zeros from train set ----
  data.train <- rsample::training(dataset)
  train.w <- data.train %>%
    dplyr::pull(w)
  train.sim <- data.train %>%
    dplyr::pull(sim)
  train.dist <- data.train %>%
    dplyr::pull(dist)
  ntrain <- data.train %>%
    nrow()
  ## Create distance in the w-dist-sim space the training set ----
  dist3d_train <- disma3d(
    train.w %>%
      standardized(),
    train.dist, train.sim, ntrain
  ) %>%
    unlist() %>%
    pracma::Reshape(ntrain, ntrain)
  hclust_train <- hclust(dist3d_train %>%
    as.dist(), method = "complete"
  )
  ## Get id from train set ----
  k <- 2
  train.ids <- cutree(hclust_train, k = k)
  ## Take out zeros from the test set ----
  data.test <- rsample::testing(dataset)
  test.set <- data.test#[zeros$test_zeros == "1", ]
  test.sim <- test.set %>%
    dplyr::pull(sim)
  test.dist <- test.set %>%
    dplyr::pull(dist)
  #*** Predict w test from sim test v
  test.w <- wlm_weight_pred(
    train.w,
    train.sim,
    train.dist,
    test.sim,
    test.dist
  )
  ntest <- test.set %>%
    nrow()
  ## Create distance matrix between train and test set ----
  dist3d.train.test <- disma3d_truncated(
    train.w %>%
      standardized(),
    train.dist,
    train.sim,
    test.w %>%
      standardized(),
    test.dist,
    test.sim,
    ntrain, ntest
  )
  ## Assign membership to test set ----
  test.ids.mat <- constrained_hclust(
    dist3d.train.test, train.ids, ntrain, ntest, k
  ) %>%
    unlist() %>%
    pracma::Reshape(ntest, k) %>%
    t()
  test.ids <- rep(0, ntest)
  for (i in 1:k)
    test.ids[test.ids.mat[i, ] == 1] <- i
  ## Create distance matrix between train and test set ----
  dist2list <- distma2centroid(
    train.dist, train.sim, test.dist, test.sim,
    train.ids, test.ids, ntrain, ntest, k
  )
  ## Put centroids relative distances into dataframes
  # dist2centroid <- dist2list[[1]] %>%
  #   unlist() %>%
  #   pracma::Reshape(max(c(ntrain, ntest)), 2 * k) %>%
  #   t()
  # k.train <- dist2centroid[1:k, ] %>%
  #   t()
  # colnames(k.train) <- paste("id", 1:k, sep = "_")
  # k.train <- k.train %>%
  #   dplyr::as_tibble()
  # k.train <- k.train[k.train$id_1 != 0, ]
  # k.test <- dist2centroid[(k + 1):(2 * k), ] %>%
  #   t()
  # colnames(k.test) <- paste("id", 1:k, sep = "_")
  # k.test <- k.test %>%
  #   dplyr::as_tibble()
  # k.test <- k.test[k.test$id_1 != 0, ]
  # # Take out possible NAN columns
  # k.train <- k.train[, !is.nan(colSums(k.train))]
  # k.test <- k.test[, !is.nan(colSums(k.test))]
  # # Assign relative distances to train & test set ----
  # data.train <- data.train %>%
  #   dplyr::bind_cols(k.train)
  # test.set <- test.set %>%
  #   dplyr::bind_cols(k.test)
  # ntest <- test.set %>%
  #   nrow()
  data.train <- data.train %>%
    dplyr::bind_cols(
      dplyr::tibble(
        id = train.ids
      )
    )
  test.set <- test.set %>%
    dplyr::bind_cols(
      dplyr::tibble(
        id = test.ids
      )
    )
  ntest <- test.set %>%
    nrow()
  ## Plot dis-sim train  communities ----
  K.coords <- dist2list[[2]] %>%
    unlist() %>%
    pracma::Reshape(2, k) %>%
    t()
  K.coords <- K.coords[!is.nan(rowSums(K.coords)), ]
  make_linkcommunities(
    data.train, test.set, K.coords, train.ids, k, serie, inst,
    subfolder = subfolder, on = T
  )
  # Save data ----
  data <- structure(
    list(data = data.train %>% dplyr::bind_rows(test.set),
         in_id = 1:ntrain,
         out_id = (ntrain + 1):(ntrain + ntest)
    ),
    class = "rsplit"
  )
  return(data)
}