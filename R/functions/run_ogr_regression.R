just_distance <- function(dataset) {
  train <- rsample::training(dataset) %>%
    dplyr::select(-sim)
  test <- rsample::testing(dataset) %>%
    dplyr::select(-sim)
  data <- structure(
    list(
      data = train %>%
        rbind(test) %>%
        dplyr::as_tibble(),
      in_id = nrow(train) %>%
        seq_len(),
      out_id = (nrow(train) + 1):(nrow(train) + nrow(test))
    ),
    class = "rsplit"
  )
  return(data)
}

run_ogr_regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_ogr.R")
    source("functions/latent_distance_process.R")
    source("functions/latent_distance_process_matrix.R")
    print("* Starting xgb regression")
    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
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
          net, nt, nodes, labels, split, inst, zero = T
        )
        ids <- datasets$ids
        st <- datasets$st
        datasets <- datasets$data
        datasets <- latent_distance_process_matrix(
          datasets, st, serie, inst,
          save = F
        )
        source("functions/plot_latent_process.R")
        plot_latent_process(
          datasets, serie, inst, subfolder = "VariableEvolution"
        )
        source("functions/plot_latent_process_b.R")
        plot_latent_process_b(
          datasets, serie, inst, subfolder = "VariableEvolution"
        )
        datasets <- datasets$data
        ### Run LM ----
        lm_model <- parsnip::linear_reg(
          mode = "regression",
          engine = "lm"
        )
        source("functions/lm_model.R")
        lm_model <- lm_model(lm_model, datasets)
        # ### Save ----
        save_work_ogr(
          lm_model, datasets, ids, inst, mats,
          serie = serie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}

set_subserie <- function(subserie, n=100, m=10) {
  long_serie <- 1:n
  return(
    long_serie[(((subserie - 1) * m) + 1):(subserie * m)]
  )
}

run_ogr_regression_crc <- function(
  net, nt, nodes, labels, serie, subserie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_ogr_crc.R")
    source("functions/latent_distance_process.R")
    source("functions/latent_distance_process_matrix.R")
    print("* Starting xgb regression")
    options(mc.cores = parallel::detectCores())
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
        datasets <- get.dataset.sim_dist(
          net, nt, nodes, labels, split, inst, zero = T
        )
        ids <- datasets$ids
        st <- datasets$st
        datasets <- datasets$data
        datasets <- latent_distance_process_matrix(
          datasets, st, serie, inst,
          save = T
        )
        datasets <- datasets$data
        ### Run LM ----
        lm_model <- parsnip::linear_reg(
          mode = "regression",
          engine = "lm"
        )
        source("functions/lm_model.R")
        lm_model <- lm_model(lm_model, datasets)
        # ### Save ----
        save_work_ogr_crc(
          lm_model, datasets, ids, inst, mats,
          serie = serie, subserie = subserie, trial = i, fold = j
        )
      }
    }
  } else
    print("** No xgb regressions")
}