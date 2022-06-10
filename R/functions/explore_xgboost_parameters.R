explore.xgboost.parameters <- function(
  model, datasets, mats, ids, inst,
  subfolder="", suffix="", on=T) {
  if (on) {
    source("functions/eval_tools.R")
    train_pred <- model$train_pred
    test_pred <- model$test_pred
    # Get predictions from train set ----
    train_prediction <- train_pred %>%
      tune::collect_predictions() %>%
      dplyr::group_by(id)
    train_prediction <- train_prediction %>%
      dplyr::select(w, .pred)
    # RMAE train set ----
    train_scores <- train_pred %>%
      tune::collect_metrics()
    train_rmae <- train_scores$mean[1]
    train_rmse <- train_scores$mean[2]
    dplyr::tibble(
      rmae = train_rmae,
      rmse = train_rmse
    ) %>%
      knitr::kable() %>%
      print()
    # Explainer train ----
    explainer <- DALEXtra::explain_xgboost(
      model = model$train_fit,
      data = rsample::training(datasets),
      y = rsample::training(datasets)$w,
      label = "XGboost"
    )
    mstudio <- modelStudio::modelStudio(explainer)
    r2d3::save_d3_html(
      mstudio,
      file = "%s/%s/Regression/XGBOOST/%s/model_train_%s.html" %>%
        sprintf(inst$plot, inst$folder, subfolder, suffix)
    )
    # Explainer test ----
    explainer <- DALEX::explain(
      model = model$train_fit,
      data = rsample::testing(datasets),
      y = rsample::testing(datasets)$w,
      label = "XGboost"
    )
    mstudio <- modelStudio::modelStudio(explainer)
    r2d3::save_d3_html(
      mstudio,
      file = "%s/%s/Regression/XGBOOST/%s/model_test_%s.html" %>%
        sprintf(inst$plot, inst$folder, subfolder, suffix)
    )
    # Check in test ----
    # Get predictions from test set ----
    test_prediction <- test_pred %>%
      tune::collect_predictions() %>%
      dplyr::select(w, .pred)
     # RMAE test set ----
    test_scores <- test_pred %>%
      tune::collect_metrics()
    test_rmse <- test_scores$.estimate[1]
    test_rmae <- test_scores$.estimate[2]
    dplyr::tibble(
      rmae = test_rmae,
      rmse = test_rmse
    ) %>%
      knitr::kable() %>%
      print()
    # Information ----
    # Graphs ----
    ## Residuals ----
    source("functions/plot_diagnosis.R")
    print("* Residuals train")
    train_residuals <- train_prediction %>%
      dplyr::arrange(.pred) %>%
      dplyr::mutate(.resid = (w - .pred)) %>%
      dplyr::mutate(.stdresid = .resid / sd(.resid))
    p <- plot.diagnosis(train_residuals)
    purb <- ggpubr::ggarrange(
      p$res.fit, p$norm.qq, p$scale.location, p$den.res,
      nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")
    )
    png(
      "%s/%s/Regression/XGBOOST/%s/residuals_train%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, suffix),
      width = 9, height = 8, res = 200, units = "in"
    )
    print(purb)
    dev.off()
    print("* Residuals test")
    test_residuals <- test_prediction %>%
      dplyr::filter(w > 0) %>%
      dplyr::arrange(.pred) %>%
      dplyr::mutate(.resid = (w - .pred)) %>%
      dplyr::mutate(.stdresid = .resid / sd(.resid))
    p <- plot.diagnosis(test_residuals)
    purb <- ggpubr::ggarrange(
      p$res.fit, p$norm.qq, p$scale.location, p$den.res,
      nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")
    )
    png(
      "%s/%s/Regression/XGBOOST/%s/residuals_test%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, suffix),
      width = 9, height = 8, res = 200, units = "in"
    )
    print(purb)
    dev.off()
    ## VIP
    p.vip <- test_pred %>%
    purrr::pluck(".workflow", 1) %>%
    tune::extract_fit_parsnip() %>%
    vip::vip()
    png(
      "%s/%s/Regression/XGBOOST/%s/vip_test%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, suffix),
      width = 6, height = 5, res = 200, units = "in"
    )
    print(p.vip)
    dev.off()
  }
}