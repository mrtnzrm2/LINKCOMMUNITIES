save.work.experiment <- function(model, datasets, ids, serie, inst, mats, suffix=""){
  source("functions/eval_tools.R")
  preprocessing.recipe <- recipes::recipe(w ~ ., rsample::training(datasets)) %>% recipes::prep()
  xgboost.model <- model$model %>%
    tune::finalize_model(model$parameters)
  # Compute in train ----
  train.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::training(datasets))
  train.model <- xgboost.model %>%
    # fit the model on all the training data
    parsnip::fit(
      formula = w ~ ., 
      data    = train.processed
    ) 
  best.score <- xgboost::xgb.importance(model = train.model$fit) %>% dplyr::as_tibble()
  train.prediction <- train.model %>%
    predict(new_data = train.processed) %>%
    dplyr::bind_cols(rsample::training(datasets))
  # Compute in test ----
  test.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::testing(datasets))
  test.prediction <- train.model %>%
    predict(new_data = test.processed) %>%
    dplyr::bind_cols(rsample::testing(datasets))
  # Ensemble parameters ----
  source("functions/rmae.R")
  parameters <- model$parameters %>% 
    dplyr::bind_cols(
      dplyr::tibble(
        train.rmae = rmae(train.prediction$w, train.prediction$.pred),
        test.rmae = rmae(test.prediction$w %>% get.nonzero(mats$test), test.prediction$.pred %>% get.nonzero(mats$test)),
        norm.trarin.rmae = rmae(train.prediction$w, norm.pred.train(train.prediction$.pred, mats$train, ids$train)),
        norm.test.rmae = rmae(test.prediction$w %>% get.nonzero(mats$test), norm.pred.test(test.prediction$.pred, mats$test, ids$test))
      ), best.score
    )
  # Ensemble regression
  regression.values <- dplyr::tibble(w=rsample::testing(datasets)$w %>% get.nonzero(mats$test), pred=test.prediction$.pred %>% get.nonzero(mats$test), norm.pred=norm.pred.test(test.prediction$.pred, mats$test, ids$test), id=ids$test, serie=paste(serie, suffix, sep = "_"))
  # Check if files exist ----
  if (!file.exists("../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>% sprintf(inst$folder, inst$common, serie)))
    write.csv(regression.values, "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>% sprintf(inst$folder, inst$common, serie), row.names = F)
  else{
    r <- read.csv("../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>% sprintf(inst$folder, inst$common, serie))
    regression.values <- r %>% dplyr::bind_rows(regression.values)
    write.csv(regression.values, "../CSV/%s/XGBOOST/%s/model_predictions_%i.csv" %>% sprintf(inst$folder, inst$common, serie), row.names = F)
  }
  if (!file.exists("../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>% sprintf(inst$folder, inst$common, serie)))
    write.csv(parameters, "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>% sprintf(inst$folder, inst$common, serie), row.names = F)
  else{
    r <- read.csv("../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>% sprintf(inst$folder, inst$common, serie))
    parameters <- r %>% dplyr::bind_rows(parameters)
    write.csv(parameters, "../CSV/%s/XGBOOST/%s/xgboost_parameters_%i.csv" %>% sprintf(inst$folder, inst$common, serie), row.names = F)
  }
}