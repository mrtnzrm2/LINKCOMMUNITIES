get.optimized.xgboost <- function(datasets, inst, grid.size=60, filename="1", save=T){
  if (save){
    # Preprocessing ----
    preprocessing.recipe <- recipes::recipe(
      w ~ ., 
      rsample::training(datasets)
    ) %>%
      recipes::prep()
    train.cv.folds <- recipes::bake(
      preprocessing.recipe, 
      new_data = rsample::training(datasets)
    ) %>% 
      rsample::vfold_cv(v = 5)
    # XGBoost model specification ----
    xgboost.model <- parsnip::boost_tree(
      mode = "regression",
      trees = 1000,
      min_n = e1071::tune(),
      tree_depth = e1071::tune(),
      learn_rate = e1071::tune(),
      loss_reduction = e1071::tune()
    ) %>% 
      parsnip::set_engine("xgboost",
      objective = "reg:squarederror",
      lambda = 20
    )
    # Grid specification ----
    xgboost.params <- dials::parameters(
        dials::min_n(),
        dials::tree_depth(),
        dials::learn_rate(),
        dials::loss_reduction()
    )
    xgboost.grid <- dials::grid_max_entropy(
        xgboost.params,
        size = grid.size
    )
    # Define the workflow ----
    xgboost.wf <- workflows::workflow() %>%
      workflows::add_model(xgboost.model) %>%
      workflows::add_formula(w ~ .)
    # Hyperparameter tuning ----
    source("functions/num-rmae.R")
    xgboost.tuned <- tune::tune_grid(
      object = xgboost.wf,
      resamples = train.cv.folds,
      grid = xgboost.grid,
      metrics = yardstick::metric_set(rmae, yardstick::rmse),
      control = tune::control_grid(verbose = TRUE)
    )
    # Review hyperparameters ----
    p <- xgboost.tuned %>%
      tune::show_best(metric = "rmae") %>%
      knitr::kable()
    print(p)
    p <- xgboost.tuned %>%
      tune::show_best(metric = "rmse") %>%
      knitr::kable()
    print(p)
    # Select best parameters ----
    xgboost.best.params <- xgboost.tuned %>%
      tune::select_best("rmse")
    # Save object ----
    l <- list(model=xgboost.model, parameters=xgboost.best.params, workflow=xgboost.wf, tune=xgboost.tuned)
    # saveRDS(l, "../RDS/imputation/%s/XGBOOST/%s.rds" %>% sprintf(inst$common, filename))
  } else
    l <- readRDS("../RDS/imputation/%s/XGBOOST/%s.rds" %>% sprintf(inst$common, filename))
  
  return(l)
}