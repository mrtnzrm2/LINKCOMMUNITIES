explore.xgboost.parameters <- function(model, datasets, mats, ids, inst, subfolder="", suffix="", on=T){
  if (on){
    source("functions/rmae.R")
    source("functions/eval_tools.R")
    # Check train set ----
    preprocessing.recipe <- recipes::recipe(w ~ ., rsample::training(datasets)) %>% recipes::prep()
    train.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::training(datasets))
    train.model <- model %>%
      parsnip::fit(
        formula = w ~ ., 
        data    = train.processed
      ) 
    train.prediction <- train.model %>%
      predict(new_data = train.processed) %>%
      dplyr::bind_cols(rsample::training(datasets))
    # RMAE train set ----
    train.score <- dplyr::tibble(rmae=rmae(train.prediction$w, train.prediction$.pred),
                                 norm.trarin.rmae = rmae(train.prediction$w, norm.pred.train(train.prediction$.pred, mats$train, ids$train))) %>%
      knitr::kable()
    print(train.score)
    # Explainer train ----
    explainer <- DALEX::explain(
      model=train.model,
      data=rsample::training(datasets),
      y= rsample::training(datasets)$w,
      label="XGboost"
    )
    mstudio <- modelStudio::modelStudio(explainer)
    r2d3::save_d3_html(mstudio, file = "%s/%s/Regression/XGBOOST/%s/model_train_%s.html" %>% sprintf(inst$plot, inst$folder, subfolder, suffix))
    # Explainer test ----
    explainer <- DALEX::explain(
      model=train.model,
      data=rsample::testing(datasets),
      y= rsample::testing(datasets)$w,
      label="XGboost"
    )
    mstudio <- modelStudio::modelStudio(explainer)
    r2d3::save_d3_html(mstudio, file = "%s/%s/Regression/XGBOOST/%s/model_test_%s.html" %>% sprintf(inst$plot, inst$folder, subfolder, suffix))
    # Check in test ----
    test.processed <- recipes::bake(preprocessing.recipe,  new_data = rsample::testing(datasets))
    test.prediction <- train.model %>%
      predict(new_data = test.processed) %>%
      dplyr::bind_cols(rsample::testing(datasets))
    
    w.test <- test.prediction$w %>% get.nonzero(mats$test)
    pred.test <- test.prediction$.pred %>% get.nonzero(mats$test)
    pred.test.norm <- test.prediction$.pred %>% norm.pred.test(mats$test, ids$tes)
    
    test.score <- dplyr::tibble( test.rmae = rmae(w.test, pred.test),
                                 norm.test.rmae = rmae(w.test, pred.test.norm)) %>%
      knitr::kable()
    print(test.score)
    # Graphs ----
    source("functions/plot_diagnosis.R")
    print("* Residuals train")
    train.residuals <- train.prediction %>%
      dplyr::arrange(.pred) %>%
      dplyr::mutate(.resid = (w - .pred)) %>%
      dplyr::mutate(.stdresid = .resid/sd(.resid))
    p <- plot.diagnosis(train.residuals)
    purb <- ggpubr::ggarrange(p$res.fit, p$norm.qq, p$scale.location, p$den.res, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
    png("%s/%s/Regression/XGBOOST/%s/residuals_train%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, suffix), width = 9, height = 8, res = 200, units = "in")
    print(purb)
    dev.off()
    
    print("* Residuals test")
    test.residuals.norm <- dplyr::tibble(
      w=w.test,
      .pred=pred.test,
      .resid=(w.test - pred.test),
      .stdresid=(w.test - pred.test)/sd(w.test - pred.test))
    p <- plot.diagnosis(test.residuals.norm)
    purb <- ggpubr::ggarrange(p$res.fit, p$norm.qq, p$scale.location, p$den.res, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
    png("%s/%s/Regression/XGBOOST/%s/residuals_test%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, suffix), width = 9, height = 8, res = 200, units = "in")
    print(purb)
    dev.off()
    
  }
}