save_work_ogr_crc <- function(
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
      "../RDS/imputation/%s/XGBOOST/%i/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
      )
    )
      saveRDS(
        lm_review,
        "../RDS/imputation/%s/XGBOOST/%i/model_parameters_%i.rds" %>%
        sprintf(inst$common, serie, subserie)
      )
  else{
    r <- readRDS(
      "../RDS/imputation/%s/XGBOOST/%i/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
    )
    lm_review <- r %>%
      dplyr::bind_rows(lm_review)
    saveRDS(
      lm_review,
      "../RDS/imputation/%s/XGBOOST/%i/model_parameters_%i.rds" %>%
      sprintf(inst$common, serie, subserie)
    )
  }
}