save_work_ogr <- function(
  model, datasets, ids, inst, mats,
  serie=0, trial=0, fold=0) {
  source("functions/eval_tools.R")
  source("functions/num-rmae.R")
  train_pred <- model$train_pred
  test_pred <- model$test_pred
  source("functions/rmae.R")
  train_rmae <- rmae(rsample::training(datasets)$w, train_pred$.pred)
  test_rmae <- rmae(rsample::testing(datasets)$w, test_pred$.pred)
  print("RMAE ---> Train: %.3f    Test: %.3f" %>%
   sprintf(
     train_rmae,
     test_rmae
   )
  )
  # Ensemble regression predictions ----
  zeros <- rsample::testing(datasets)$w == 0
  regression_values <- dplyr::tibble(
    w = rsample::testing(datasets)$w[!zeros],
    pred = test_pred$.pred[!zeros],
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
      "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie)
      )
    )
      saveRDS(
        regression_values,
        "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
        sprintf(inst$common, serie)
      )
  else{
    r <- readRDS(
      "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie)
    )
    regression_values <- r %>%
      dplyr::bind_rows(regression_values)
    saveRDS(
      regression_values,
      "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>%
      sprintf(inst$common, serie)
    )
  }
  ## For linear models, get fit parameters and statistics
  lm_review <- model$train_fit %>%
    broom::tidy() %>%
    dplyr::bind_cols(
      dplyr::tibble(
        serie = serie,
        trial = trial,
        fold = fold
      )
    )
  if (!file.exists(
      "../RDS/imputation/%s/XGBOOST/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie)
      )
    )
      saveRDS(
        lm_review,
        "../RDS/imputation/%s/XGBOOST/model_parameters_%i.rds" %>%
        sprintf(inst$common, serie)
      )
  else{
    r <- readRDS(
      "../RDS/imputation/%s/XGBOOST/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie)
    )
    lm_review <- r %>%
      dplyr::bind_rows(lm_review)
    saveRDS(
      lm_review,
      "../RDS/imputation/%s/XGBOOST/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie)
    )
  }
}