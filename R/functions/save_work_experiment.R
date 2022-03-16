save.work.experiment <- function(
  model, datasets, ids, inst, mats,
  serie=0, trial=0, fold=0) {
  source("functions/eval_tools.R")
  train_pred <- model$train_pred
  test_pred <- model$test_pred
  best_score <- xgboost::xgb.importance(
            model = train_pred$fit) %>%
            dplyr::as_tibble()
  test_prediction <- test_pred %>%
    tune::collect_predictions() %>%
    dplyr::group_by(id)
  train_rmae <- train_pred %>%
    tune::collect_metrics()
  train_rmae <- train_rmae$mean[1]
  # Ensemble parameters ----
  source("functions/rmae.R")
  parameters <- model$parameters %>%
    dplyr::bind_cols(
      dplyr::tibble(
        train.rmae = train_rmae,
        rmae(rsample::testing(datasets)$w, test_prediction$.pred),
      ), best_score
    )
  # Ensemble regression
  regression_values <- dplyr::tibble(
    w = rsample::testing(datasets)$w %>%
     get.nonzero(mats$test),
    pred = test.prediction$.pred %>%
      get.nonzero(mats$test),
    id = ids$test, serie = serie, trial = trial, fold = fold
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
      regression.values,
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