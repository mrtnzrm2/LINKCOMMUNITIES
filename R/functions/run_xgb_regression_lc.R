run_xgb_regression_lc <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_regression_parsnip.R")
    Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
    Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
    print("* Starting xgb regression")
    # Start experiment ----
    for (i in 1:100) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        print("** Trial: %i kfold: %i" %>%
          sprintf(i, j)
        )
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
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
        ### Set-up XGBOOST ----
        xgboost_model <- 
        # parsnip::linear_reg(penalty = 100) %>%
        #   parsnip::set_engine("glmnet")
        parsnip::boost_tree(
          mode = "regression",
          trees = 500,
          min_n = 26,
          tree_depth =  10,
          learn_rate = 0.0043587,
          loss_reduction = 0.0007983
          ) %>%
            parsnip::set_engine(
              "xgboost",
              objective = "reg:squarederror",
              lambda = 0
            )
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ### Plot RSME
        make_train_rmse(
          xgboost_model$train_fit, inst,
          subfolder = "Residuals",
          suffix = "_trial_%i_fold_%i" %>% sprintf(i, j),
          on = F
        )
        ### Explore model ----
        explore.xgboost.parameters(
          xgboost_model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>% as.character(), on = F
        )
        ### Save ----
        save_work_regression_parsnip(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}