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

plot.xgb.residuals <- function(serie, inst, subfolder="", suffix=""){
  source("functions/rmae.R")
  dir.create("%s/%s/Regression/XGBOOST/%s/%i" %>% sprintf(inst$plot, inst$folder, subfolder, serie), showWarnings = F)
  # Load data ----
  data <- readRDS("../RDS/imputation/%s/XGBOOST/model_predictions_%i.rds" %>% sprintf(inst$common, serie))
  # Marginalize over ids
  data <- data %>% dplyr::group_by(id) %>% dplyr::summarise(w=unique(w), stdist=unique(dist), .pred=mean(pred))
  data <- data %>% dplyr::mutate(.resid=(w-.pred), .stdresid = (w-.pred)/sd((w-.pred)), .rmae=abs(w-.pred)/w)
  source("functions/plot_diagnosis.R")
  # Plot diagnosis ----
  p <- plot.diagnosis(data)
  purb <- ggpubr::ggarrange(p$res.fit, p$norm.qq, p$scale.location, p$den.res, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
  png("%s/%s/Regression/XGBOOST/%s/%i/residuals_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 9, height = 8, res = 200, units = "in")
  print(purb)
  dev.off()
  # Plot RMAE ----
  p.wrmae <- ggplot2::ggplot(data, ggplot2::aes(w, .rmae))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("w")+
    ggplot2::ylab("Relative prediction error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange")+
    ggplot2::theme_bw()
  p.distrmae <- ggplot2::ggplot(data, ggplot2::aes(stdist, .rmae))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Standardize distance")+
    ggplot2::ylab("Relative prediction error")+
    ggplot2::stat_smooth(method = "loess", color="orange", fill="orange")+
    ggplot2::theme_bw()
  purb <- ggpubr::ggarrange(p.wrmae, p.distrmae, p.rmse, nrow = 1, ncol = 2, labels = c("A", "B" , "C"))
  png("%s/%s/Regression/XGBOOST/%s/%i/rmae_marg_id%s.png" %>% sprintf(inst$plot, inst$folder, subfolder, serie, suffix), width = 8, height = 3.5, res = 200, units = "in")
  print(purb)
  dev.off()
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
  # Run ----
  run.xgb.regression(net, nt, nodes, labels, serie, inst, kfold=3, work=T)
  # Plots ----
  make.experiment.parameters(serie, inst, subfolder="SerieAnalysis", on=F)
  # plot.xgb.residuals(serie, inst, subfolder="Residuals")
  
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
