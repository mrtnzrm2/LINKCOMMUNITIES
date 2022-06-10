multilevel_model <- "
data {
  // global
  int<lower=1> K;
  int<lower=1> groups;
  // train
  int<lower=1> ntrn;
  int<lower=1,upper=groups> group_train[ntrn];
  matrix[ntrn, K] x_train;
  vector[ntrn] y_train;
  // test
  int<lower=1> ntst;
  int<lower=1,upper=groups> group_test[ntst];
  matrix[ntst, K] x_test;
  vector[ntst] y_test;
}
parameters {
  matrix[K, groups] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0, upper=pi() / 2>[K] tau_unif;
  real<lower=0> sigma;
}
transformed parameters {
  vector<lower=0>[K] tau = 2.5 * tan(tau_unif);
  matrix[K, groups] beta = diag_pre_multiply(tau, L_Omega) * z;
}
model {
  vector[ntrn] mu_train;
  for(i in 1:ntrn) {
    mu_train[i] = x_train[i, ] * beta[, group_train[i]];
  }
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  y_train ~ normal(mu_train, sigma);
}
generated quantities {
  real prediction_train[ntrn];
  real prediction_test[ntst];
  vector[ntrn] rmae_train = rep_vector(0, ntrn);
  vector[ntst] rmae_test = rep_vector(0, ntst);
  real train_rmae = 0;
  real test_rmae = 0;
  int nonzero_test = 0;
  for (i in 1:ntrn) {
    prediction_train[i] = normal_rng(
      x_train[i,] * beta[, group_train[i]], sigma
    );
    if (y_train[i] > 0)
      rmae_train[i] = fabs(y_train[i] - prediction_train[i]) / y_train[i];
  }
  for (i in 1:ntst) {
    prediction_test[i] = normal_rng(
      x_test[i,] * beta[, group_test[i]], sigma
    );
    if (y_test[i] > 0)
      rmae_test[i] = fabs(y_test[i] - prediction_test[i]) / y_test[i];
  }
  train_rmae = mean(rmae_train);
  for (i in 1:ntst) {
    if (rmae_test[i] > 0) {
      test_rmae += rmae_test[i];
      nonzero_test += 1;
    }
  }
  if (nonzero_test > 0)
    test_rmae /= nonzero_test;
}"

stan_regression <- function(data) {
  # Separate data ----
  train_data <- rsample::training(data)
  test_data <- rsample::testing(data)
  ntrn <- train_data %>%
    nrow()
  ntst <- test_data %>%
    nrow()
  train_w <- train_data %>%
    dplyr::pull(w)
  train_sim <- train_data %>%
    dplyr::pull(sim)
  train_dist <- train_data %>%
    dplyr::pull(dist)
  train_group <- train_data %>%
    dplyr::pull(id)
  test_w <- test_data %>%
    dplyr::pull(w)
  test_sim <- test_data %>%
    dplyr::pull(sim)
  test_dist <- test_data %>%
    dplyr::pull(dist)
  test_group <- test_data %>%
    dplyr::pull(id)
  groups <- train_group %>%
    unique() %>%
    length()
  # Format data ----
  x_train <- matrix(
    c(
      rep(1, ntrn),
      train_sim,
      train_dist
    ),
    ncol = 3
  )
  x_test <- matrix(
    c(
      rep(1, ntst),
      test_sim,
      test_dist
    ),
    ncol = 3
  )
  stan_set <- list(
    K = 3,
    groups = groups,
    ntrn = ntrn,
    ntst = ntst,
    group_train = train_group,
    group_test = test_group,
    x_train = x_train,
    x_test = x_test,
    y_train = train_w,
    y_test = test_w
  )
  options(mc.cores = parallel::detectCores())
  stan_sample <- rstan::stan(
    model_code = multilevel_model,
    model_name = "multilevel_model",
    data = stan_set,
    iter = 200,
    verbose = T
  )
  print(stan_sample, pars = c("train_rmae", "test_rmae", "beta"))
}

run_hierarchical_regression <- function(
  net, nt, nodes, labels, serie, inst,
  kfold=3, work=T) {
  if (work) {
    source("functions/split_kfold.R")
    source("functions/setup_datasets.R")
    source("functions/explore_xgboost_parameters.R")
    source("functions/link_zeros.R")
    source("functions/link_classification.R")
    source("functions/save_work_regression.R")
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
        stan_regression(datasets)
      }
    }
  } else
    print("** No xgb regressions")
}