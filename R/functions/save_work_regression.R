save.work.regression <- function(model, datasets, ids, inst, mats, serie=0, trial=0, fold=0){
  source("functions/eval_tools.R")
  preprocessing.recipe <- recipes::recipe(w ~ ., rsample::training(datasets)) %>% recipes::prep()
  # Compute in train ----
  train.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::training(datasets))
  train.model <- model %>%
    # fit the model on all the training data
    parsnip::fit(
      formula = w ~ ., 
      data    = train.processed
    ) 
  train.prediction <- train.model %>%
    predict(new_data = train.processed) %>%
    dplyr::bind_cols(rsample::training(datasets))
  # Compute in test ----
  test.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::testing(datasets))
  test.prediction <- train.model %>%
    predict(new_data = test.processed) %>%
    dplyr::bind_cols(rsample::testing(datasets))
  # Print RMAE ----
  source("functions/rmae.R")
  train.rmae <- rmae(train.prediction$w, train.prediction$.pred)
  test.rmae <- rmae(test.prediction$w %>% get.nonzero(mats$test), test.prediction$.pred %>% get.nonzero(mats$test))
  print("**** RMAE ---> train: %.3f test: %.3f" %>% sprintf(train.rmae, test.rmae))
  # Ensemble regression
  regression.values <-
    dplyr::tibble(
      w=rsample::testing(datasets)$w %>%
        get.nonzero(mats$test),
      pred=test.prediction$.pred %>% 
        get.nonzero(mats$test),
      dist=rsample::testing(datasets)$dist %>%
        get.nonzero(mats$test),
      sim=rsample::testing(datasets)$sim %>%
        get.nonzero(mats$test),
      norm.pred=norm.pred.test(test.prediction$.pred, mats$test, ids$test),
      id=ids$test,
      serie=serie,
      trial=trial,
      fold=fold
      )
  # Check if files exist ----
  if (!file.exists("../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie)))
    saveRDS(regression.values, "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie))
  else{
    r <- readRDS("../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie))
    regression.values <- r %>% dplyr::bind_rows(regression.values)
    saveRDS(regression.values, "../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie))
  }
}