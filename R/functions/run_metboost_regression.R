run_metboost_regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/link_classification.R")
    source("functions/save_work_regression_metboost.R")
    Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
    Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
    print("* Starting metboost regression")
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
        ### Link classification ----
        datasets <- link.classification(
          datasets, "test_zeros", serie, inst,
          subfolder = "VariableEvolution"
        )
        # Start model ----
        print("** Run metb")
        model <- mvtboost::metb(
          y = rsample::training(datasets)$w,
          X = rsample::training(datasets) %>%
            dplyr::select(-w),
          id = "id",
          n.trees = 200,
          cv.folds = 5,
          n.minobsinnode = 26,
          shrinkage = 0.0043587,
          interaction.depth = 10,
          mc.cores = 4,
          save.mods = T
        )
        ### Save ----
        save_work_regression_metboost(
          model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}