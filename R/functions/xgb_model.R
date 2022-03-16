xgb_model <- function(model, datasets) {
    source("functions/num-rmae.R")
    # For performance ----
    ## CV partition ----
    cv_folds <- rsample::vfold_cv(
        rsample::training(datasets),
        v = 10,
        strata = w
    )
    ## Prepare recipe ----
    model_recipe <- recipes::recipe(w ~ ., datasets)
    ## Compute in train ----
    train_wf <- workflows::workflow() %>%
        workflows::add_model(model) %>%
        workflows::add_recipe(model_recipe)
    train_pred <- train_wf  %>%
        tune::fit_resamples(
            resamples = cv_folds,
            metrics = yardstick::metric_set(
                yardstick::rmse,
                rmae
            ),
        control = tune::control_resamples(save_pred = T)
        )
    ## Compute in test ----
    test_pred <- tune::last_fit(
        train_wf,
        split = datasets,
        metrics = yardstick::metric_set(
            yardstick::rmse,
            rmae
        )
    )
    # Fit ----
    train_model <- model %>%
      parsnip::fit(
        formula = w ~ .,
        data    = rsample::training(datasets)
      )
    return(
        list(
            train_pred = train_pred,
            test_pred = test_pred,
            train_fit = train_model
        )
    )
}