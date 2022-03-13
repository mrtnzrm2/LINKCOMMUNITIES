make.train.rmse <- function(datasets, model, inst, subfolder="", suffix="", on=T){
  if (on){
    print("*** Plot train rmse")
    # Plot importance ----
    preprocessing.recipe <- recipes::recipe(w ~ ., rsample::training(datasets))
    final.xgboost <- tune::last_fit(model, preprocessing.recipe, datasets)
    p <- final.xgboost %>% tune::extract_fit_parsnip()
    p <- p$fit
    p.dat <- ggplot2::ggplot(p$evaluation_log, ggplot2::aes(iter, training_rmse))+
      ggplot2::geom_line()
    
    png("%s/%s/Regression/XGBOOST/%s/eval_niter%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, suffix), width = 6, height = 5, res = 200, units = "in")
    print(p.dat)
    dev.off()
  } else
    print("*** No train rmse")
}

make.tree.depth <- function(tuned, inst, subfolder="", filename="", on=T){
  if (on){
    print("*** Plotting tree_depth")
    p <- tuned %>%
      tune::collect_metrics() %>%
      dplyr::mutate(tree_depth = factor(tree_depth)) %>%
      ggplot2::ggplot(ggplot2::aes(min_n, mean, color = tree_depth)) +
      ggplot2::geom_line(size = 0.5, alpha = 0.6) +
      ggplot2::geom_point(size = 1) +
      ggplot2::facet_wrap(~ .metric, scales = "free", nrow = 3)
    png("%s/%s/Regression/XGBOOST/%s/tree_depth_%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, filename), width = 6, height = 6, res = 200, units = "in")
    print(p)
    dev.off()
  } else
    print("*** No tree_depth")
}

make.experiment.parameters <- function(series, inst, subfolder="", on=T){
  if (on){
    source("functions/plot_experiment_performance.R")
    print("* Plotting XGBOOST parameters")
    dir.create(sprintf('%s/%s/%s', inst$plot, inst$mainfolder, subfolder), showWarnings = F)
    plot.experiment.performance(series, inst, subfolder=subfolder)
  } else
    print("* No XGBOOST parameters")
}

make.xgb.residuals <- function(series, inst, subfolder="", on=T){
  if (on){
    source("functions/plot_xgb_residuals.R")
    print("** Plotting xgb residuals")
    plot.xgb.residuals(serie, inst, subfolder=subfolder)
  } else
    print("** No xgb residuals")
}

make.parameters <- function(features, inst, subfolder="", on=T){
  if (on){
    print("*** Printing parameters plots")
    source('functions/plot_process_parameters.R')
    plot.process.parameters(inst$plot, inst$folder, features, subfolder=subfolder)  
  } else
    print("*** No parameters")
}

run.xgb.exp.regression <- function(net, nt, nodes, labels, serie, inst, kfold=3, work=T){
  if (work){
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/get_optimized_xgboost.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_experiment.R")
    print("** Starting experiment")
    # Start experiment ----
    for (i in 1:100){
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold){
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ## Get datasets ----
        datasets <- get.dataset(net, nt, nodes, labels, split, inst, filename=serie %>% as.character(), save=T)
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        ## Hyperparameter optimization and Model ----
        xmodel <- get.optimized.xgboost(datasets, inst, grid.size = 60, filename =serie %>% as.character(), save=T)
        xgboost.model <- xmodel$model %>% tune::finalize_model(xmodel$parameters)
        ## Plots ----
        make.train.rmse(datasets, xgboost.model, inst, subfolder="Residuals_EXP", on=T)
        explore.xgboost.parameters(xgboost.model, datasets, mats, ids, inst, subfolder="Residuals_EXP", suffix=serie %>% as.character(), on=F)
        ## Save ----
        save.work.experiment(xmodel, datasets, ids, serie, inst, mats, suffix = i)
      }
    }
  } else
    print("** No experiment")
}

make.linkcommunities <- function(train, test, serie, inst, subfolder="", on=T){
  if (on){
    print("*** Plotting link communities in the dis-sim space")
    source("functions/plot_linkcomm_validation.R")
    plot.linkcommm.validation(train, test, serie, inst, subfolder=subfolder)
  } else
    print("*** No dis-sim space")
}

run.xgb.regression <- function(net, nt, nodes, labels, serie, inst, kfold=3, work=T){
  if (work){
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression.R")
    print("* Starting xgb regression")
    # Start experiment ----
    for (i in 1:100){
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:kfold){
        print("** Trial: %i kfold: %i" %>% sprintf(i,j))
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset(net, nt, nodes, labels, split, inst, filename=serie %>% as.character(), save=T)
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
        ### Plots ----
        make.train.rmse(datasets, xgboost.model, inst, subfolder="Residuals", suffix="_trial_%i_fold_%i" %>% sprintf(i,j), on=F)
        explore.xgboost.parameters(xgboost.model, datasets, mats, ids, inst, subfolder="Residuals", suffix=serie %>% as.character(), on=F)
        # ### Save ----
        save.work.regression(xgboost.model, datasets, ids, inst, mats, serie=serie, trial=i, fold=j)
      }
    }
  } else
    print("** No xgb regressions")
}

link.classification <- function(dataset, serie, inst, subfolder=""){
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>% sprintf(inst$plot, inst$folder, subfolder, serie), showWarnings = F)
  source("functions/df_to_adj.R")
  train.w <- rsample::training(dataset) %>% dplyr::pull(w)
  train.sim <- rsample::training(dataset)  %>% dplyr::pull(sim)
  train.dist <- rsample::training(dataset)  %>% dplyr::pull(dist)
  ntrain <- rsample::training(dataset) %>% nrow()
  test.sim <- rsample::testing(dataset)  %>% dplyr::pull(sim)
  test.dist <- rsample::testing(dataset)  %>% dplyr::pull(dist)
  ntest <- rsample::testing(dataset) %>% nrow()
  # Create distance from the training set ----
  Rcpp::sourceCpp("../cpp/distance_matrix.cpp")
  dist3d.train <- disma3d(train.w, train.dist, train.sim, ntrain) %>% unlist() %>% pracma::Reshape(ntrain, ntrain)
  hclust.train <- hclust(dist3d.train %>% as.dist() , method = "complete")
  ## Process hclust ----
  source('functions/process_hclust_distance.R')
  hclust.process <- process.hclust.dist(hclust.train)
  ## Plot SS parameter ----
  make.parameters(hclust.process, inst,  subfolder=subfolder, on=T)
  # Assign id to train set ----
  k <- hclust.process %>% dplyr::filter(SS == max(SS, na.rm = T)) %>% dplyr::pull(K) %>% unique() %>% sort()
  if (length(k) > 1)
    k <- k[1]
  data.train <- rsample::training(dataset) %>% dplyr::bind_cols(dplyr::tibble(id=cutree(hclust.train, k=k)))
  ## Create distance matrix between train and test set ----
  dist2d.train.test <- disma2d_truncated(train.dist, train.sim, test.dist, test.sim, ntrain, ntest)
  ## Find final cluster ----
  Rcpp::sourceCpp("../cpp/constrained_hcluster.cpp")
  test.communities <- constrained_hclust(dist2d.train.test, data.train$id, ntrain, ntest, k) %>% unlist() %>% pracma::Reshape(ntest, k) %>% t()
  # Assign id to test set ----
  data.test <- rsample::testing(dataset)
  data.test$id <- 0
  for (i in 1:k)
    data.test$id[test.communities[i,] == 1] <- i
  data.train$id <- data.train$id %>% as.character()
  data.test$id <- data.test$id %>% as.character()
  ## Plot dis-sim train-test communities ----
  make.linkcommunities(data.train, data.test, serie, inst, subfolder = subfolder, on=T)
  # Save data ----
  data <- structure(
    list(data = data.train %>% dplyr::bind_rows(data.test),
         in_id = 1:ntrain,
         out_id = (ntrain+1):(ntrain+ntest)
    ),
    class = "rsplit"
  )
  data.train$cat <- "train"
  data.test$cat <- "test"
  dat <- data.train %>% dplyr::bind_rows(data.test)
  p <- ggplot2::ggplot(dat, ggplot2::aes(w))+
    ggplot2::facet_grid(cat ~ id)+
    ggplot2::geom_histogram()
  png("%s/%s/Regression/XGBOOST/%s/%i/hist_id.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie), width = 15, height =10 , res = 200, units = "in")
  print(p)
  dev.off()
  return(data)
}

run.xgb.regression.lc <- function(net, nt, nodes, labels, serie, inst, kfold=3, work=T){
  if (work){
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_regression.R")
    print("* Starting xgb regression")
    # Start experiment ----
    for (i in 1:1){
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:1){
        print("** Trial: %i kfold: %i" %>% sprintf(i,j))
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ### Get datasets ----
        datasets <- get.dataset(net, nt, nodes, labels, split, inst, filename=serie %>% as.character(), save=T)
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        datasets <- link.classification(datasets, serie, inst, subfolder = "VariableEvolution")
        ### Set-up XGBOOST ----
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
        # ### Save ----
        save.work.regression(xgboost.model, datasets, ids, inst, mats, serie=serie, trial=i, fold=j)
      }
    }
  } else
    print("** No xgb regressions")
}

run.xgb.exp.regression.lc <- function(net, nt, nodes, labels, serie, inst, kfold=3, work=T){
  if (work){
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/get_optimized_xgboost.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/save_work_experiment.R")
    print("** Starting experiment")
    # Start experiment ----
    for (i in 1:1){
      ## Kfold loop ----
      split.kfold <- split.dataset.kfold(nt, kfold = kfold)
      for (j in 1:1){
        ### Assemble split ----
        split <- assemble.split(split.kfold, j)
        ## Get datasets ----
        datasets <- get.dataset(net, nt, nodes, labels, split, inst, filename=serie %>% as.character(), save=T)
        ids <- datasets$ids
        mats <- datasets$matids
        datasets <- datasets$data
        datasets <- link.classification(datasets, serie, inst, subfolder = "VariableEvolution")
        ## Hyperparameter optimization and Model ----
        xmodel <- get.optimized.xgboost(datasets, inst, grid.size = 10, filename =serie %>% as.character(), save=T)
        xgboost.model <- xmodel$model %>% tune::finalize_model(xmodel$parameters)
        ## Plots ----
        make.train.rmse(datasets, xgboost.model, inst, subfolder="Residuals_EXP", on=F)
        explore.xgboost.parameters(xgboost.model, datasets, mats, ids, inst, subfolder="Residuals_EXP", suffix=serie %>% as.character(), on=T)
        ## Save ----
        save.work.experiment(xmodel, datasets, ids, serie, inst, mats, suffix = i)
      }
    }
  } else
    print("** No experiment")
}

main <- function(inst){
  library(magrittr)
  set.seed(NULL)
  
  # Get series number ----
  args = commandArgs(trailingOnly=TRUE)
  serie <- args[1] %>% as.numeric()
  
  # Parallelization ----
  all.cores <- parallel::detectCores(logical = FALSE)
  doParallel::registerDoParallel(cores = all.cores-1)
  
  nlog10 <- T
  nt <- 50   ## Number of  columns
  nodes <- 107
  
  # Load network ----
  source('functions/load_net.R')
  netx <- load.net(inst)
  net <- netx$net
  leaves <- netx$leaves
  labels <- netx$labels
  source("functions/format_labels.R")
  labels <- labels %>% format.labels()
  
  # Create folder ----
  inst$mainfolder <- paste(inst$folder, "Regression", sep = "/")
  dir.create(sprintf('%s/%s', inst$plot, inst$mainfolder), showWarnings = F)
  
  if (nlog10){
    net$weight[net$weight != 0] <- log10(net$weight[net$weight != 0]) + 7
  }
  
  # Experiment ----
  run.xgb.exp.regression(net, nt, nodes, labels, serie, inst, work=F)
  run.xgb.exp.regression.lc(net, nt, nodes, labels, serie, inst, work=T)
  # Run ----
  run.xgb.regression(net, nt, nodes, labels, serie, inst, kfold=3, work=F)
  run.xgb.regression.lc(net, nt, nodes, labels, serie, inst, kfold=3, work=F)
  # Plots ----
  make.experiment.parameters(serie, inst, subfolder="SerieAnalysis", on=F)
  make.xgb.residuals(serie, inst, subfolder = "Residuals", on=F)
  
}


### HEAD ###
folder <- 'merged'
distances <- 'original'
model <- 'normal'
csv.name <- 'fln'
labels.name <- 'rlabels'  # imputation_labels  rlabels

common.path <- paste(distances, model, sep = '/')
csv.path <- paste(folder, 'imputation', common.path, paste0(csv.name,'.csv'), sep = '/') # merged/imputation/tracto2016/zz_model/
labels.path <- paste(folder, 'labels', common.path, paste0(labels.name,'.csv'), sep = '/')
plot.path <- '../plots'

path.list <- list(csv=csv.path, 
                  labels=labels.path, 
                  plot=plot.path, 
                  common=common.path,
                  folder=folder,
                  model=model,
                  distances=distances)

main(path.list)
