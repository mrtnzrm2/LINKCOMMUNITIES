save_work_regression_metboost <- function(
  model, datasets, ids, inst, mats,
  serie=0, trial=0, fold=0) {
  source("functions/eval_tools.R")
  source("functions/num-rmae.R")
  train_pred <- model %>%
    predict(
      rsample::training(datasets) %>%
        dplyr::select(-w),
      "id"
    )
  test_pred <- model %>%
    predict(
      rsample::testing(datasets) %>%
        dplyr::select(-w),
      "id"
    )
  source("functions/rmae.R")
  test_rmae <- rmae(rsample::testing(datasets)$w, test_pred$yhat)
  train_rmae <- rmae(rsample::training(datasets)$w, train_pred$yhat)
  print("RMAE ---> Train: %.3f    Test: %.3f" %>%
   sprintf(
     train_rmae,
     test_rmae
   )
  )
  # # Ensemble regression
  # zeros <- rsample::testing(datasets)$w == 0
  # regression_values <- dplyr::tibble(
  #   w = rsample::testing(datasets)$w[!zeros],
  #   pred = test_prediction$.pred[!zeros],
  #   dist = rsample::testing(datasets)$dist[!zeros],
  #   sim = rsample::testing(datasets)$sim[!zeros],
  #   train_rmae = train_rmae,
  #   test_rmae = test_rmae,
  #   id = ids$test[!zeros],
  #   serie = serie,
  #   trial = trial,
  #   fold = fold
  # )
  # # Check if files exist ----
  # if (!file.exists(
  #     "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
  #     sprintf(inst$common, serie)
  #     )
  #   )
  #     saveRDS(
  #       regression_values,
  #       "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
  #       sprintf(inst$common, serie)
  #     )
  # else{
  #   r <- readRDS(
  #     "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
  #     sprintf(inst$common, serie)
  #   )
  #   regression_values <- r %>%
  #     dplyr::bind_rows(regression_values)
  #   saveRDS(
  #     regression_values,
  #     "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
  #     sprintf(inst$common, serie)
  #   )
  # }
}