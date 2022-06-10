data {
  //*** others v
  int<lower=1> d;
  int<lower=2> K;
  int<lower=1> M;
  //*** train v
  int<lower=1> ntrn;
  int<lower=1, upper=K> category_train[ntrn];
  real dist_train[ntrn];
  real sim_train[ntrn];
  int source_train[ntrn];
  int target_train[ntrn];
  //*** test v
  int<lower=1> ntst;
  real dist_test[ntst];
  real sim_test[ntst];
  int source_test[ntst];
  int target_test[ntst];
}

parameters {
  real alpha;
  real<lower=machine_precision()> sigma;
  real<lower=machine_precision()> rho;
  real<lower=machine_precision()> rho_s;
  ordered[K-1] b;
  vector<lower=-1, upper=1>[d] z[M];
  vector<lower=-1, upper=1>[d] z_s[M];
}

transformed parameters {
  vector<lower=0>[ntrn + ntst] asym;
  //** Estimating categories for train set v
  for(i in 1:ntrn) {
    asym[i] = sqrt(
      dot_self(z[target_train[i]] - z[source_train[i]]) +
      dot_self(z[target_train[i]] - z_s[source_train[i]])
    );
  }
  //** Estimating categories for test set v
  for(i in 1:ntst) {
    asym[i + ntrn] = sqrt(
      dot_self(z[target_test[i]] - z[source_test[i]]) +
      dot_self(z[target_test[i]] - z_s[source_test[i]])
    );
  }
}

model {
  for (i in 1:M) {
    z_s[i] ~ normal(0,rho_s);
    z[i] ~ normal(0, rho);
  }
  {
    vector[K] f[ntrn];
    //** Estimating categories for train set v
    for(i in 1:ntrn) {
      f[i, 1] = 1 - Phi(
        (dist_train[i] + alpha * sim_train[i] + asym[i] - b[1]) / sigma
      );
      for(j in 2:(K-1))
        f[i, j] = Phi(
          (dist_train[i] + alpha * sim_train[i] + asym[i] - b[j-1]) / sigma
        ) -
          Phi(
            (dist_train[i] + alpha * sim_train[i] + asym[i] - b[j]) / sigma
          );
      f[i, K] = Phi(
        (dist_train[i] + alpha * sim_train[i] + asym[i] - b[K-1] ) / sigma
      );
      category_train[i] ~ categorical(f[i]);
    }
  }
}

generated quantities{
  vector[ntrn] train_ldist;
  vector[ntst] test_ldist;
  for (i in 1:ntrn)
    train_ldist[i] = (dist_train[i] + alpha * sim_train[i] + asym[i]) / sigma;
  for (i in 1:ntst)
    test_ldist[i] = (dist_test[i] + alpha * sim_test[i] + asym[i + ntrn]) / sigma;
}
