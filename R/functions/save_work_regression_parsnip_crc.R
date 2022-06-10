save_work_regression_parsnip_crc <- function(
  model, datasets, ids, inst, mats,
  serie=0, subserie=0, trial=0, fold=0) {
  if (!file.exists("../RDS/imputation/%s/XGBOOST/%i" %>%
        sprintf(inst$common, serie)
      )
    )
    dir.create(
      "../RDS/imputation/%s/XGBOOST/%i" %>%
        sprintf(inst$common, serie)
    )
  source("functions/eval_tools.R")
  source("functions/num-rmae.R")
  train_pred <- model$train_pred
  test_pred <- model$test_pred
  test_prediction <- test_pred %>%
    tune::collect_predictions() %>%
    dplyr::group_by(id)
  train_rmae <- train_pred %>%
    tune::collect_metrics()
  train_rmae <- train_rmae$mean[1]
  source("functions/rmae.R")
  test_rmae <- rmae(rsample::testing(datasets)$w, test_prediction$.pred)
  print("RMAE ---> Train: %.3f    Test: %.3f" %>%
   sprintf(
     train_rmae,
     test_rmae
   )
  )
  # Ensemble regression
  zeros <- rsample::testing(datasets)$w == 0
  regression_values <- dplyr::tibble(
    w = rsample::testing(datasets)$w[!zeros],
    pred = test_prediction$.pred[!zeros],
    dist = rsample::testing(datasets)$dist[!zeros],
    # sim = rsample::testing(datasets)$sim[!zeros],
    train_rmae = train_rmae,
    test_rmae = test_rmae,
    id = ids$test[!zeros],
    serie = serie,
    trial = trial,
    fold = fold
  )
  # Check if files exist ----
  if (!file.exists(
      "../RDS/imputation/%s/XGBOOST/%i/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
      )
    )
      saveRDS(
        regression_values,
        "../RDS/imputation/%s/XGBOOST/%i/model_predictions_%i.rds" %>%
        sprintf(inst$common, serie, subserie)
      )
  else{
    r <- readRDS(
      "../RDS/imputation/%s/XGBOOST/%i/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
    )
    regression_values <- r %>%
      dplyr::bind_rows(regression_values)
    saveRDS(
      regression_values,
      "../RDS/imputation/%s/XGBOOST/%i/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
    )
  }
}