save.work.experiment <- function(
  model, datasets, ids, inst, mats,
  serie=0, trial=0, fold=0) {
  source("functions/eval_tools.R")
  train_pred <- model$train_pred
  test_pred <- model$test_pred
  train_fit <- model$train_fit
  best_score <- xgboost::xgb.importance(
            model = train_fit$fit) %>%
            dplyr::as_tibble()
  test_prediction <- test_pred %>%
    tune::collect_predictions() %>%
    dplyr::group_by(id)
  train_rmae <- train_pred %>%
    tune::collect_metrics()
  train_rmae <- train_rmae$mean[1]
  # Ensemble parameters ----
  source("functions/rmae.R")
  test_rmae <- rmae(rsample::testing(datasets)$w, test_prediction$.pred)
  parameters <- model$parameters %>%
    dplyr::bind_cols(
      dplyr::tibble(
        train.rmae = train_rmae,
        test.rmae = test_rmae,
      ), best_score
    )
  # Ensemble regression
  zeros <- rsample::testing(datasets)$w == 0
  regression_values <- dplyr::tibble(
    w = rsample::testing(datasets)$w[!zeros],
    pred = test_prediction$.pred[!zeros],
    dist = rsample::testing(datasets)$dist[!zeros],
    sim = rsample::testing(datasets)$sim[!zeros],
    train_rmae = train_rmae,
    test_rmae = test_rmae,
    id = ids$test[!zeros],
    serie = serie,
    trial = trial,
    fold = fold
  )
  # Check if files exist ----
  if (!file.exists(
        "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>%
          sprintf(inst$folder, inst$common, serie)
      )
  )
    write.csv(
      regression_values,
      "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie), row.names = F
    )
  else{
    r <- read.csv(
      "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie)
    )
    regression_values <- r %>%
      dplyr::bind_rows(regression_values)
    write.csv(
      regression_values,
      "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie), row.names = F
    )
  }
  if (!file.exists(
        "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>%
          sprintf(inst$folder, inst$common, serie)
      )
  )
    write.csv(
      parameters,
      "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie), row.names = F
    )
  else{
    r <- read.csv(
      "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie)
    )
    parameters <- r %>%
      dplyr::bind_rows(parameters)
    write.csv(
      parameters,
      "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>%
        sprintf(inst$folder, inst$common, serie), row.names = F
    )
  }
}