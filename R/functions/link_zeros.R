link.zero <- function(dataset, serie, inst, subfolder="") {
  dir.create(
        "%s/%s/Regression/XGBOOST/%s/%i" %>%
            sprintf(inst$plot, inst$folder, subfolder, serie),
        showWarnings = F
    )
  data.train <- rsample::training(dataset)
  data.test <- rsample::testing(dataset)
  ## Get important variables ----
  train.w <- data.train %>%
    dplyr::pull(w)
  train.sim <- data.train %>%
    dplyr::pull(sim)
  train.dist <- data.train %>%
    dplyr::pull(dist)
  ntrain <- data.train %>%
    nrow()
  test.sim <- data.test %>%
    dplyr::pull(sim)
  test.dist <- data.test %>%
    dplyr::pull(dist)
  ntest <- data.test %>%
    nrow()
  ## Create distance in the w-dist-sim space the training set ----
  dist3d.train <- disma3d(
    ifelse(train.w > 0, 3, 3),
    train.dist, train.sim, ntrain
  ) %>%
    unlist() %>%
    pracma::Reshape(ntrain, ntrain)
  hclust.train <- hclust(
    dist3d.train %>%
      as.dist(),
    method = "ward.D"
  )
  ## Get id from train set ----
  k <- 2
  train.ids <- cutree(hclust.train, k = k)
  ## Create distance matrix between train and test set ----
  dist2list <- distma2zero(
    train.dist, train.sim, test.dist, test.sim,
    train.ids, ntrain, ntest, k
  )
  ## Put centroids relative distances into dataframes ----
  dist2centroid <- dist2list[[1]] %>%
    unlist() %>%
    pracma::Reshape(max(c(ntrain, ntest)), 2 * k) %>%
    t()
  k.train <- dist2centroid[1:k, ] %>%
    t() %>%
    dplyr::as_tibble()
  colnames(k.train) <- paste("id", 1:k, sep = "_")
  k.train <- k.train[k.train$id_1 != 0, ]
  k.test <- dist2centroid[(k + 1):(2 * k), ] %>%
    t() %>%
    dplyr::as_tibble()
  colnames(k.test) <- paste("id", 1:k, sep = "_")
  ### Take out rows in exces ----
  k.test <- k.test[k.test$id_1 != 0, ]
  ## Tak position of the 1 and 0 center of mass ----
  K_coords <- dist2list[[2]] %>%
    unlist() %>%
    pracma::Reshape(2, k) %>%
    t()
  ## Set-up XGBOOST ----
  ### Prepare train dataset for classification ----
  train.bin <- data.train
  train.bin$w  <- ifelse(train.bin$w > 0, "1", "0")
  train.bin$w <- factor(train.bin$w, levels = c("1", "0"))
  ### Prepare test dataset for classification ----
  test.bin <- data.test
  test.bin$w  <- ifelse(test.bin$w > 0, "1", "0")
  test.bin$w <- factor(test.bin$w, levels = c("1", "0"))
  ### Put centroid distances ----
  train.bin <- train.bin %>%
    dplyr::bind_cols(k.train)
  test.bin <- test.bin  %>%
    dplyr::bind_cols(k.test)
  ### Merge both sets in a unique split
  data <- structure(
    list(data = train.bin %>%
      dplyr::bind_rows(test.bin),
      in_id = 1:ntrain,
      out_id = (ntrain + 1):(ntrain + ntest)
    ),
    class = "rsplit"
  )
  ### CV validation ----
  cv_folds <- rsample::vfold_cv(
    train.bin,
    v = 10,
    strata = w
  )
  ### Creating recipe ----
  xgboost.recipe <- recipes::recipe(w ~ ., train.bin)
  ### Create the xgboost ----
  train.model <- parsnip::boost_tree(
      mode = "classification",
      trees = 1000,
      min_n = 50,
      tree_depth =  7,
      learn_rate = 0.01,
      loss_reduction = 0.01
  ) %>%
    parsnip::set_engine("xgboost", objective = "binary:logistic")
  ### Create workflow ----
  xgboost.wf <- workflows::workflow() %>%
    workflows::add_recipe(xgboost.recipe) %>%
    workflows::add_model(train.model)
  # Compute in train ----
  train.model <-  xgboost.wf %>%
    tune::fit_resamples(
      resamples = cv_folds,
      metrics = yardstick::metric_set(
        yardstick::roc_auc,
        yardstick::accuracy,
        yardstick::precision),
      control = tune::control_resamples(save_pred = T)
  )
  ### Evaluation of the train set ----
  print(train.model %>%
    tune::collect_metrics(summarize = T)
  )
  ### Collect train set predictions ----
  train.prediction <- train.model %>%
    tune::collect_predictions()
  print(train.prediction %>%
    yardstick::conf_mat(w, .pred_class)
  )
  ### ROC ----
  p.roc <- train.prediction %>%
    dplyr::group_by(id) %>%
    yardstick::roc_curve(w, .pred_1) %>%
    ggplot2::autoplot()
  ### Density graphs ----
  p.den <- train.prediction %>%
    ggplot2::ggplot(ggplot2::aes(x = .pred_1, fill = w)) +
    ggplot2::geom_density(alpha = 0.5)
  ### Confusion matrix ----
  p.conf <- train.model %>%
    tune::collect_predictions() %>%
    yardstick::conf_mat(w, .pred_class) %>%
    ggplot2::autoplot(type = "heatmap")
  ### Plot
  p <- cowplot::plot_grid(p.roc, p.den, p.conf, ncol = 3)
  png("%s/%s/Regression/XGBOOST/%s/%i/binary_train_eval_id.png"
    %>% sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 13, height = 5, res = 200, units = "in")
  print(p)
  dev.off()
  # Compute in test set ----
  test.model <- tune::last_fit(xgboost.wf,
    split = data,
    metrics = yardstick::metric_set(
      yardstick::accuracy,
      yardstick::roc_auc,
      yardstick::precision)
  )
  ### Evaluation of the test set ----
  print(test.model %>%
    tune::collect_metrics()
  )
  ### Test set plots
  ### Vip
  p.vip <- test.model %>%
    purrr::pluck(".workflow", 1) %>%
    tune::extract_fit_parsnip() %>%
    vip::vip()
  ### ROC
  p.roc <- test.model %>%
    tune::collect_predictions() %>%
    yardstick::roc_curve(w, .pred_1) %>%
    ggplot2::autoplot()
  ### Confusion matrix
  p.conf <- test.model %>%
    tune::collect_predictions() %>%
    yardstick::conf_mat(w, .pred_class) %>%
    ggplot2::autoplot(type = "heatmap")
  ### Plot
  p <- cowplot::plot_grid(p.vip, p.roc,p.conf, ncol = 3)
  png("%s/%s/Regression/XGBOOST/%s/%i/binary_test_eval_id.png"
    %>% sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 12, height = 5, res = 200, units = "in")
  print(p)
  dev.off()
  # Save
  df <- list(
    test_zeros <- test.model %>%
      tune::collect_predictions() %>%
      dplyr::group_by(id) %>%
      dplyr::pull(.pred_class),
    K_coords = K_coords
  )
  return(df)
}