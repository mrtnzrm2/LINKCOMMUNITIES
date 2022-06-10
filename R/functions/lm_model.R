lm_model <- function(model, datasets) {
    source("functions/num-rmae.R")
    # For performance ----
    train_set <- rsample::training(datasets)
    print(train_set)
    train_model <- model %>%
      parsnip::fit(
        formula = w ~ .,
        data    = train_set
      )
    return(
        list(
            train_pred = predict(
              train_model,
              train_set
            ),
            test_pred = predict(
              train_model,
              rsample::testing(datasets)
            ),
            train_fit = train_model
        )
    )
}