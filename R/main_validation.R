make_linkcommunities <- function(
  train, test, K_coords, train_ids, K, serie, inst,
  subfolder="", on=T) {
  if (on) {
    print("*** Plotting link communities in the dis-sim space")
    source("functions/plot_linkcomm_validation.R")
    train$id <- train_ids
    plot.linkcommm.validation(train, test, K_coords, K, serie, inst,
      subfolder = subfolder
    )
  } else
    print("*** No dis-sim space")
}

make_train_rmse <- function(
  model, inst,
  subfolder="", suffix="", on=T) {
  if (on) {
    print("*** Plot train rmse")
    # Plot importance ----
    p <- model %>%
      tune::extract_fit_engine()
    p.dat <- ggplot2::ggplot(
      p$evaluation_log,
      ggplot2::aes(iter, training_rmse)
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw()
    png(
      "%s/%s/Regression/XGBOOST/%s/eval_niter_train_%s.png" %>%
      sprintf(inst$plot, inst$folder, subfolder, suffix),
      width = 6, height = 5, res = 200, units = "in"
    )
    print(p.dat)
    dev.off()
  } else
    print("*** No train rmse")
}

make_tree_depth <- function(tuned, inst, subfolder="", filename="", on=T) {
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

make_experiment_parameters <- function(series, inst, subfolder="", on=T) {
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

make_xgb_residuals <- function(series, inst, subfolder="", on=T) {
  if (on) {
    source("functions/plot_xgb_residuals.R")
    print("** Plotting xgb residuals")
    plot.xgb.residuals(series, inst, subfolder = subfolder)
  } else
    print("** No xgb residuals")
}

make_parameters <- function(features, inst, subfolder="", on=T) {
  if (on) {
    print("*** Printing parameters plots")
    source("functions/plot_process_parameters.R")
    plot.process.parameters(
      inst$plot, inst$folder, features,
      subfolder = subfolder
    )  
  } else
    print("*** No parameters")
}

make_zeros <- function(train, test, kcoords, serie, inst, subfolder="", on=T) {
  if (on) {
    print("*** Plotting zeros information")
    source("functions/plot_zeros.R")
    # Some plots ----
    plot.zeros(train, test, kcoords, serie, inst, subfolder = subfolder)
  } else
    print("*** No zeros information")
}

r2Dcpp <- function(A, nc) {
  cxx.mat <- A[, 1] %>%
    list()
  for (i in 2:nc)
    cxx.mat <- c(cxx.mat, A[, i] %>% list())
  return(cxx.mat)
}

main <- function(inst) {
  library(magrittr)
  # set.seed(12345)
  # setwd("/Users/jmarti53/Documents/ND/LINKCOMMUNITIES/MAC/220126/R")
  # Get series number ----
  args <- commandArgs(trailingOnly = TRUE)
  serie <- args[1] %>%
    as.numeric()
  subserie <- args[2] %>%
    as.numeric()
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
  source("functions/run_xgb_exp_regression.R")
  run_xgb_exp_regression(
    net, nt, nodes, labels, serie, inst,
    work = F
  )
  source("functions/run_xgb_exp_regression_lc.R")
  run_xgb_exp_regression_lc(
    net, nt, nodes, labels, serie, inst,
    work = F
  )
  # Run ----
  source("functions/run_xgb_regression.R")
  run_xgb_regression(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = F
  )
  source("functions/run_xgb_regression_lc.R")
  run_xgb_regression_lc(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = F
  )
  source("functions/run_hierarchical_regression.R")
  run_hierarchical_regression(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = F
  )
  source("functions/run_metboost_regression.R")
  run_metboost_regression(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = F
  )
  source("functions/run_ogr_regression.R")
  run_ogr_regression(
    net, nt, nodes, labels, serie, inst,
    kfold = 3, work = T
  )
  run_ogr_regression_crc(
    net, nt, nodes, labels, serie, subserie, inst,
    kfold = 5, work = F
  )
  # Plots ----
  make_experiment_parameters(
    serie, inst,
    subfolder = "SerieAnalysis", on = F
  )
  make_xgb_residuals(
    serie, inst,
    subfolder = "Residuals", on = F
  )
}

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
