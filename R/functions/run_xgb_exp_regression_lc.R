run_xgb_exp_regression_lc <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/get_optimized_xgboost.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_experiment.R")
    Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
    Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
    print("** Starting experiment")
    # Start experiment ----
    for (i in 1:1) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:1) {
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ## Get datasets ----
        datasets <- get.dataset.sim_dist(
          net, nt, nodes, labels, split, inst
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        # test_zeros <- link.zero(
        #   datasets, serie, inst,
        #   subfolder = "VariableEvolution"
        # )
        datasets <- link.classification(
          datasets, "test_zeros", serie, inst,
          subfolder = "VariableEvolution"
        )
        ## Hyperparameter optimization and Model ----
        xmodel <- get.optimized.xgboost(
          datasets, inst,
          grid.size = 60, filename = serie %>% as.character(), save = T
        )
        xgboost_model <- xmodel$model %>%
          tune::finalize_model(xmodel$parameters)
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ## Plots ----
        make_train_rmse(
          xgboost_model$train_fit, inst,
          subfolder = "Residuals_EXP", on = F
        )
        explore.xgboost.parameters(
          xgboost_model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>% as.character(), on = F
        )
        ## Save ----
        save.work.experiment(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No experiment")
}