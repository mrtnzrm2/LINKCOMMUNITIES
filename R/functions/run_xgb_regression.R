set_subserie <- function(subserie, n=100, m=10) {
  long_serie <- 1:n
  return(
    long_serie[(((subserie - 1) * m) + 1):(subserie * m)]
  )
}

run_xgb_regression_crc <- function(
  net, nt, nodes, labels, serie, subserie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression_parsnip_crc.R")
    print("* Starting xgb regression")
    # Start experiment ----
    for (i in set_subserie(subserie)) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        print("** Trial: %i kfold: %i" %>%
          sprintf(i, j)
        )
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset.fln_dist(
          net, nt, nodes, labels, split, inst
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        ### Run XGBOOST ----
        xgboost_model <- parsnip::boost_tree(
            mode = "regression",
            trees = 1000, # 750 1000
            min_n = 20, # 18 20
            tree_depth =  15, # 2 15
            learn_rate = 0.0142631, # 0.0094047 0.0142631
            loss_reduction = 0.1355005 # 0.0458003 0.1355005
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        # ### Save ----
        save_work_regression_parsnip_crc(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, subserie = subserie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}

run_xgb_regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression_parsnip.R")
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
        datasets <- get.dataset.fln_dist(
          net, nt, nodes, labels, split, inst
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        ### Run XGBOOST ----
        xgboost_model <- parsnip::boost_tree(
            mode = "regression",
            trees = 1000, # 750 1000
            min_n = 20, # 18 20
            tree_depth =  15, # 2 15
            learn_rate = 0.0142631, # 0.0094047 0.0142631
            loss_reduction = 0.1355005 # 0.0458003 0.1355005
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ### Plots ----
        make_train_rmse(
          xgboost_model$train_fit, inst,
          subfolder = "Residuals_EXP",
          suffix = "_trial_%i_fold_%i" %>% sprintf(i, j),
          on = F
        )
        explore.xgboost.parameters(
          xgboost_model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>% as.character(), on = F
        )
        # ### Save ----
        save_work_regression_parsnip(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}