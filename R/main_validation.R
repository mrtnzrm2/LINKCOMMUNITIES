make.train.rmse <- function(
  model, inst,
  subfolder="", suffix="", on=T) {
  if (on) {
    print("*** Plot train rmse")
    # Plot importance ----
    p <- model$test %>%
      tune::extract_fit_parsnip()
    p <- p$fit
    p.dat <- ggplot2::ggplot(
      p$evaluation_log,
      ggplot2::aes(iter, training_rmse)
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw()
    png(
      "%s/%s/Regression/XGBOOST/%s/eval_niter%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, suffix),
      width = 6, height = 5, res = 200, units = "in"
    )
    print(p.dat)
    dev.off()
  } else
    print("*** No train rmse")
}

make.tree.depth <- function(tuned, inst, subfolder="", filename="", on=T) {
  if (on) {
    print("*** Plotting tree_depth")
    p <- tuned %>%
      tune::collect_metrics() %>%
      dplyr::mutate(tree_depth = factor(tree_depth)) %>%
      ggplot2::ggplot(ggplot2::aes(min_n, mean, color = tree_depth)) +
      ggplot2::geom_line(size = 0.5, alpha = 0.6) +
      ggplot2::geom_point(size = 1) +
      ggplot2::facet_wrap(~ .metric, scales = "free", nrow = 3)
    png(
      "%s/%s/Regression/XGBOOST/%s/tree_depth_%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, filename),
      width = 6, height = 6, res = 200, units = "in"
    )
    print(p)
    dev.off()
  } else
    print("*** No tree_depth")
}

make.experiment.parameters <- function(series, inst, subfolder="", on=T) {
  if (on) {
    source("functions/plot_experiment_performance.R")
    print("* Plotting XGBOOST parameters")
    dir.create(
      sprintf(
        "%s/%s/%s", inst$plot, inst$mainfolder, subfolder
        ),
      showWarnings = F
    )
    plot.experiment.performance(series, inst, subfolder = subfolder)
  } else
    print("* No XGBOOST parameters")
}

make.xgb.residuals <- function(series, inst, subfolder="", on=T) {
  if (on){
    source("functions/plot_xgb_residuals.R")
    print("** Plotting xgb residuals")
    plot.xgb.residuals(series, inst, subfolder = subfolder)
  } else
    print("** No xgb residuals")
}

make.parameters <- function(features, inst, subfolder="", on=T) {
  if (on){
    print("*** Printing parameters plots")
    source("functions/plot_process_parameters.R")
    plot.process.parameters(
      inst$plot, inst$folder, features,
      subfolder = subfolder
    )  
  } else
    print("*** No parameters")
}

make.linkcommunities <- function(
  train, test, K.coords, train.ids, K, serie, inst,
  subfolder="", on=T) {
  if (on ) {
    print("*** Plotting link communities in the dis-sim space")
    source("functions/plot_linkcomm_validation.R")
    train$id <- train.ids
    plot.linkcommm.validation(train, test, K.coords, K, serie, inst,
      subfolder = subfolder
    )
  } else
    print("*** No dis-sim space")
}

make.zeros <- function(test, serie, inst, subfolder="", on=T) {
  if (on) {
    print("*** Plotting seros in sim-dist space")
    # Some plots ----
    test$zero <- ifelse(test$w == 0, "B", "A")
    p <- ggplot2::ggplot(test, ggplot2::aes(dist, sim, color = zero)) +
      ggplot2::geom_point(size = 0.5)
    png("%s/%s/Regression/XGBOOST/%s/%i/zerot_test.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, serie),
      width = 6, height = 5, res = 200, units = "in")
    print(p)
    dev.off()
  } else
    print("*** No dis-sim space")
}

run.xgb.exp.regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/get_optimized_xgboost.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_experiment.R")
    print("** Starting experiment")
    # Start experiment ----
    for (i in 1:100) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ## Get datasets ----
        datasets <- get.dataset(
          net, nt, nodes, labels, split, inst,
          filename = serie %>%
            as.character(),
          save = T
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        ## Hyperparameter optimization and Model ----
        xmodel <- get.optimized.xgboost(
          datasets, inst,
          grid.size = 60, filename = serie %>%
            as.character(),
          save = T
        )
        xgboost.model <- xmodel$model %>%
          tune::finalize_model(xmodel$parameters)
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ## Plots ----
        make.train.rmse(
          xgboost.model, inst,
          subfolder = "Residuals_EXP", on = T
        )
        explore.xgboost.parameters(
          xgboost.model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>%
            as.character(),
          on = F
        )
        ## Save ----
        save.work.experiment(
          xmodel, datasets, ids, serie, inst, mats,
          suffix = i
        )
      }
    }
  } else
    print("** No experiment")
}

run.xgb.regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T){
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression.R")
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
        datasets <- get.dataset(
          net,
          nt,
          nodes,
          labels,
          split,
          inst,
          filename = serie %>%
           as.character(),
          save = T
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        ### Run XGBOOST ----
        xgboost.model <-
          parsnip::boost_tree(
            mode = "regression",
            trees = 500,
            min_n = 15,
            tree_depth =  3,
            learn_rate = 0.01,
            loss_reduction = 0.01
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ### Plots ----
        make.train.rmse(
          xgboost.model, inst,
          subfolder = "Residuals",
          suffix = "_trial_%i_fold_%i" %>% sprintf(i,j),
          on = F
        )
        explore.xgboost.parameters(
          xgboost.model, datasets, mats, ids, inst,
          subfolder = "Residuals", suffix = serie %>% as.character(), on = F
        )
        # ### Save ----
        save.work.regression(
          xgboost.model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}

r2Dcpp <- function(A, nc) {
  cxx.mat <- A[, 1] %>%
    list()
  for (i in 2:nc)
    cxx.mat <- c(cxx.mat, A[, i] %>% list())
  return(cxx.mat)
}

link.zero <- function(dataset, serie, inst, subfolder="") {
  data.train <- rsample::training(dataset)
  data.test <- rsample::testing(dataset)
  ## Get important variables ----
  train.w <- data.train %>%
    dplyr::pull(w)
  train.sim <- data.train %>%
    dplyr::pull(sim)
  train.dist <- data.train %>%
    dplyr::pull(dist)
  ntrain <- data.train %>%
    nrow()
  test.sim <- data.test %>%
    dplyr::pull(sim)
  test.dist <- data.test %>%
    dplyr::pull(dist)
  ntest <- data.test %>%
    nrow()
  ## Create distance in the w-dist-sim space the training set ----
  dist3d.train <- disma3d(
    ifelse(train.w > 0, 3, 3),
    train.dist, train.sim, ntrain
  ) %>%
    unlist() %>%
    pracma::Reshape(ntrain, ntrain)
  hclust.train <- hclust(
    dist3d.train %>%
      as.dist(),
    method = "ward.D"
  )
  ## Get id from train set ----
  k <- 2
  train.ids <- cutree(hclust.train, k = k)
  ## Create distance matrix between train and test set ----
  dist2list <- distma2zero(
    train.dist, train.sim, test.dist, test.sim, train.ids, ntrain, ntest, k
  )
  ## Put centroids relative distances into dataframes
  dist2centroid <- dist2list[[1]] %>%
    unlist() %>%
    pracma::Reshape(max(c(ntrain, ntest)), 2 * k) %>%
    t()
  k.train <- dist2centroid[1:k, ] %>% t() %>% dplyr::as_tibble()
  colnames(k.train) <- paste("id", 1:k, sep = "_")
  k.train <- k.train[k.train$id_1 != 0, ]
  k.test <- dist2centroid[(k + 1):(2 * k), ] %>% t() %>% dplyr::as_tibble()
  colnames(k.test) <- paste("id", 1:k, sep = "_")
  k.test <- k.test[k.test$id_1 != 0, ]
  ## Set-up XGBOOST ----
  ### Prepare train dataset for classification ----
  train.bin <- data.train
  train.bin$w  <- ifelse(train.bin$w > 0, "1", "0")
  train.bin$w <- factor(train.bin$w, levels = c("1", "0"))
  ### Prepare test dataset for classification ----
  test.bin <- data.test
  test.bin$w  <- ifelse(test.bin$w > 0, "1", "0")
  test.bin$w <- factor(test.bin$w, levels = c("1", "0"))
  ### Put centroid distances ----
  train.bin <- train.bin %>%
    dplyr::bind_cols(k.train)
  test.bin <- test.bin  %>%
    dplyr::bind_cols(k.test)
  ### Merge both sets in a unique split
  data <- structure(
    list(data = train.bin %>%
      dplyr::bind_rows(test.bin),
      in_id = 1:ntrain,
      out_id = (ntrain + 1):(ntrain + ntest)
    ),
    class = "rsplit"
  )
  ### CV validation ----
  cv_folds <- rsample::vfold_cv(
    train.bin,
    v = 10,
    strata = w
  )
  ### Creating recipe ----
  xgboost.recipe <- recipes::recipe(w ~ ., train.bin)
  ### Create the xgboost ----
  train.model <- parsnip::boost_tree(
      mode = "classification",
      trees = 1000,
      min_n = 100,
      tree_depth =  7,
      learn_rate = 0.01,
      loss_reduction = 0.01
  ) %>%
    parsnip::set_engine("xgboost", objective = "binary:logistic")
  ### Create workflow ----
  xgboost.wf <- workflows::workflow() %>%
    workflows::add_recipe(xgboost.recipe) %>%
    workflows::add_model(train.model)
  # Compute in train ----
  train.model <-  xgboost.wf %>%
    tune::fit_resamples(
      resamples = cv_folds,
      metrics = yardstick::metric_set(
        yardstick::roc_auc,
        yardstick::accuracy,
        yardstick::precision),
      control = tune::control_resamples(save_pred = T)
  )
  ### Evaluation of the train set ----
  print(train.model %>%
    tune::collect_metrics(summarize = T)
  )
  ### Collect train set predictions ----
  train.prediction <- train.model %>%
    tune::collect_predictions()
  print(train.prediction %>%
    yardstick::conf_mat(w, .pred_class)
  )
  ### ROC ----
  p.roc <- train.prediction %>%
    dplyr::group_by(id) %>%
    yardstick::roc_curve(w, .pred_1) %>%
    ggplot2::autoplot()
  ### Density graphs ----
  p.den <- train.prediction %>%
    ggplot2::ggplot(ggplot2::aes(x = .pred_1, fill = w)) +
    ggplot2::geom_density(alpha = 0.5)
  ### Confusion matrix ----
  p.conf <- train.model %>%
    tune::collect_predictions() %>%
    yardstick::conf_mat(w, .pred_class) %>%
    ggplot2::autoplot(type = "heatmap")
  ### Plot
  p <- cowplot::plot_grid(p.roc, p.den, p.conf, ncol = 3)
  png("%s/%s/Regression/XGBOOST/%s/%i/train_eval_id.png"
    %>% sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 13, height = 5, res = 200, units = "in")
  print(p)
  dev.off()
  # Compute in test set ----
  test.model <- tune::last_fit(xgboost.wf,
    split=data,
    metrics = yardstick::metric_set(
      yardstick::accuracy,
      yardstick::roc_auc,
      yardstick::precision)
  )
  ### Evaluation of the test set ----
  print(test.model %>%
    tune::collect_metrics()
  )
  ### Test set plots
  ### Vip
  p.vip <- test.model %>%
    purrr::pluck(".workflow", 1) %>%
    tune::extract_fit_parsnip() %>%
    vip::vip()
  ### ROC
  p.roc <- test.model %>%
    tune::collect_predictions() %>%
    yardstick::roc_curve(w, .pred_1) %>%
    ggplot2::autoplot()
  ### Confusion matrix
  p.conf <- test.model %>%
    tune::collect_predictions() %>%
    yardstick::conf_mat(w, .pred_class) %>%
    ggplot2::autoplot(type = "heatmap")
  ### Plot
  p <- cowplot::plot_grid(p.vip, p.roc,p.conf, ncol = 3)
  png("%s/%s/Regression/XGBOOST/%s/%i/test_eval_id.png"
    %>% sprintf(inst$plot, inst$folder, subfolder, serie),
    width = 12, height = 5, res = 200, units = "in")
  print(p)
  dev.off()

  return(test.model %>%
    tune::collect_predictions() %>%
    dplyr::group_by(id) %>%
    dplyr::pull(.pred_class)
  )
}

link.classification <- function(dataset, tzeros, serie, inst, subfolder="") {
  source("functions/df_to_adj.R")
  # Take out zeros from train set ----
  data.train <- rsample::training(dataset)
  data.train <- data.train[data.train$w > 0, ]
  train.w <- data.train %>%
    dplyr::pull(w)
  train.sim <- data.train %>%
    dplyr::pull(sim)
  train.dist <- data.train %>%
    dplyr::pull(dist)
  ntrain <- data.train %>%
    nrow()
  ## Create distance in the w-dist-sim space the training set ----
  dist3d.train <- disma3d(
    train.w, train.dist, train.sim, ntrain
  ) %>%
    unlist() %>%
    pracma::Reshape(ntrain, ntrain)
  hclust.train <- hclust(dist3d.train %>%
    as.dist(), method = "complete"
  )
  ## Get id from train set ----
  k <- 3
  train.ids <- cutree(hclust.train, k = k)
  print("** Number of communities: %i" %>%
    sprintf(unique(train.ids) %>% length())
  )
  ## Take out zeros from the test set ----
  data.test <- rsample::testing(dataset)
  test.set <- data.test#[tzeros == "1", ]
  test.sim <- test.set %>%
    dplyr::pull(sim)
  test.dist <- test.set %>%
    dplyr::pull(dist)
  ntest <- test.set %>%
    nrow()
  ## Create distance matrix between train and test set ----
  dist2d.train.test <- disma2d_truncated(
    train.dist, train.sim, test.dist, test.sim, ntrain, ntest
  )
  ## Assign membership to test set ----
  test.ids.mat <- constrained_hclust(
    dist2d.train.test, train.ids, ntrain, ntest, k
  ) %>%
    unlist() %>%
    pracma::Reshape(ntest, k) %>%
    t()
  test.ids <- rep(0, ntest)
  for (i in 1:k)
    test.ids[test.ids.mat[i, ] == 1] <- i
  ## Create distance matrix between train and test set ----
  dist2list <- distma2centroid(
    train.dist, train.sim, test.dist, test.sim,
    train.ids, test.ids, ntrain, ntest, k
  )
  ## Put centroids relative distances into dataframes
  dist2centroid <- dist2list[[1]] %>%
    unlist() %>%
    pracma::Reshape(max(c(ntrain, ntest)), 2 * k) %>%
    t()
  k.train <- dist2centroid[1:k, ] %>%
    t()
  colnames(k.train) <- paste("id", 1:k, sep = "_")
  k.train <- k.train %>%
    dplyr::as_tibble()
  k.train <- k.train[k.train$id_1 != 0, ]
  k.test <- dist2centroid[(k + 1):(2 * k), ] %>%
    t()
  colnames(k.test) <- paste("id", 1:k, sep = "_")
  k.test <- k.test %>%
    dplyr::as_tibble()
  k.test <- k.test[k.test$id_1 != 0, ]
  # Take out possible NAN columns
  k.train <- k.train[, !is.nan(colSums(k.train))]
  k.test <- k.test[, !is.nan(colSums(k.test))]
  # Assign relative distances to train & test set ----
  data.train <- data.train %>%
    dplyr::bind_cols(k.train)
  test.set <- test.set %>%
    dplyr::bind_cols(k.test)
  ntest <- test.set %>%
    nrow()
  ## Plot dis-sim train  communities ----
  K.coords <- dist2list[[2]] %>%
    unlist() %>%
    pracma::Reshape(2, k) %>%
    t()
  K.coords <- K.coords[!is.nan(rowSums(K.coords)), ]
  make.linkcommunities(
    data.train, test.set, K.coords, train.ids, k, serie, inst,
    subfolder = subfolder, on = T
  )
  make.zeros(test.set, serie, inst, subfolder = subfolder, on = F)
  # Save data ----
  data <- structure(
    list(data = data.train %>% dplyr::bind_rows(test.set),
         in_id = 1:ntrain,
         out_id = (ntrain + 1):(ntrain + ntest)
    ),
    class = "rsplit"
  )
  return(data)
}

run_xgb_regression_lc <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression.R")
    Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
    Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
    print("* Starting xgb regression")
    # Start experiment ----
    for (i in 1:1) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:1) {
        print("** Trial: %i kfold: %i" %>%
          sprintf(i,j)
        )
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset(
          net, nt, nodes, labels, split, inst,
          filename = serie %>% as.character(), save = T
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        test_zeros <- link.zero(
          datasets, serie, inst,
          subfolder = "VariableEvolution"
        )
        datasets <- link.classification(
          datasets, test_zeros, serie, inst,
          subfolder = "VariableEvolution"
        )
        ### Set-up XGBOOST ----
        xgboost_model <-
          parsnip::boost_tree(
            mode = "regression",
            trees = 1000,
            min_n = 50,
            tree_depth =  7,
            learn_rate = 0.01,
            loss_reduction = 0.01
          ) %>%
          parsnip::set_engine("xgboost", objective = "reg:squarederror")
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ### Explore model ----
        explore.xgboost.parameters(
          xgboost_model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>% as.character(), on = T
        )
        ### Save ----
        save.work.regression(
          xgboost_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}

run_xgb_exp_regression_lc <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/get_optimized_xgboost.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_experiment.R")
    Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
    Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
    print("** Starting experiment")
    # Start experiment ----
    for (i in 1:100) {
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold) {
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ## Get datasets ----
        datasets <- get.dataset(
          net, nt, nodes, labels, split, inst,
          filename = serie %>% as.character(), save = T
        )
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        test_zeros <- link.zero(
          datasets, serie, inst,
          subfolder = "VariableEvolution"
        )
        datasets <- link.classification(
          datasets, test_zeros, serie, inst,
          subfolder = "VariableEvolution"
        )
        ## Hyperparameter optimization and Model ----
        xmodel <- get.optimized.xgboost(
          datasets, inst,
          grid.size = 5, filename = serie %>% as.character(), save = T
        )
        xgboost_model <- xmodel$model %>%
          tune::finalize_model(xmodel$parameters)
        source("functions/xgb_model.R")
        xgboost_model <- xgb_model(xgboost_model, datasets)
        ## Plots ----
        make.train.rmse(
          xgboost_model, inst,
          subfolder = "Residuals_EXP", on = T
        )
        explore.xgboost.parameters(
          xgboost_model, datasets, mats, ids, inst,
          subfolder = "Residuals_EXP", suffix = serie %>% as.character(), on = T
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

main <- function(inst) {
  library(magrittr)
  set.seed(NULL)
  # Get series number ----
  args <- commandArgs(trailingOnly = TRUE)
  # serie <- args[1] %>%
  #   as.numeric()
  serie <- 2
  # Parallelization ----
  all.cores <- parallel::detectCores(logical = FALSE)
  doParallel::registerDoParallel(cores = all.cores - 1)
  nlog10 <- T
  nt <- 50   ## Number of  columns
  nodes <- 107
  # Load network ----
  source("functions/load_net.R")
  netx <- load.net(inst)
  net <- netx$net
  labels <- netx$labels
  source("functions/format_labels.R")
  labels <- labels %>%
    format.labels()
  # Create folder ----
  inst$mainfolder <- paste(inst$folder, "Regression", sep = "/")
  dir.create(sprintf("%s/%s", inst$plot, inst$mainfolder), showWarnings = F)
  if (nlog10) {
    net$weight[net$weight != 0] <- log10(net$weight[net$weight != 0]) + 7
  }
  # Experiment ----
  run.xgb.exp.regression(
    net, nt, nodes, labels, serie, inst,
    work = F
  )
  run_xgb_exp_regression_lc(
    net, nt, nodes, labels, serie, inst,
    work = F
  )
  # Run ----
  run.xgb.regression(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = F
  )
  run_xgb_regression_lc(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = T
  )
  # Plots ----
  make.experiment.parameters(
    serie, inst,
    subfolder = "SerieAnalysis", on = F
  )
  make.xgb.residuals(
    serie, inst,
    subfolder = "Residuals", on = F
  )
}

setwd("/Users/jmarti53/Documents/ND/LINKCOMMUNITIES/MAC/220126/R")
### HEAD ###
folder <- "merged"
distances <- "original"
model <- "normal"
csv_name <- "fln"
labels_name <- "rlabels"  # imputation_labels  rlabels
common_path <- paste(distances, model, sep = "/")
# merged/imputation/tracto2016/zz_model/
csv_path <- paste(
  folder, "imputation", common_path, paste0(csv_name, ".csv"),
  sep = "/"
)
labels_path <- paste(
  folder, "labels", common_path, paste0(labels_name, ".csv"),
  sep = "/"
)
plot_path <- "../plots"
path_list <- list(csv = csv_path,
  labels = labels_path,
  plot = plot_path,
  common = common_path,
  folder = folder,
  model = model,
  distances = distances
)
main(path_list)
