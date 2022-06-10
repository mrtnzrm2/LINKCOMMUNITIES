functions {
   real latent_distance(
    real zt_x, real zt_y, real zt_z,
    real zs_x, real zs_y, real zs_z,
    real s_x, real s_y, real s_z) {
    real dist;
    dist = sqrt(
      square(zt_x - zs_x) + square(zt_y - zs_y) + square(zt_z - zs_z)
    ) +
    sqrt(
      square(zt_x - s_x) + square(zt_y - s_y) + square(zt_z - s_z)
    );
    return dist;
  }
}

data {
  int<lower=1> M;
  int<lower=1> N;
  int<lower=2> K;
  array[M,N] int<lower=1, upper=K> mat;
  array[M,N] int<lower=0, upper=1> w;
  array[M, N] real<lower=0> D;
  int<lower=1> d;
}

parameters {
    ordered[K-1] b;                                                                    
    real<lower=machine_precision()> sigma;                                            
    real<lower=machine_precision()> rho;
    real<lower=machine_precision()> rho_s;
    array[M] vector<lower=-1, upper=1>[d] z;
    array[M] vector<lower=-1, upper=1>[d] z_s;
}

transformed parameters{
  matrix<lower=0>[M, N] asym;
  // EC quadrant
  for(i in 1:N) {
    for(j in 1:N) {
      if (i < j){		  
        asym[i,j] = sqrt(dot_self(z[j]-z[i])) + sqrt(dot_self(z[j]-z_s[i]));
        asym[j,i] = sqrt(dot_self(z[i]-z[j])) + sqrt(dot_self(z[i]-z_s[j]));
      }
      if (i == j)
        asym[i, i] = 0;
    }
  }
  // Third quadrant
  for(i in (N+1):M) {
    for(j in 1:N) { 
      asym[i,j] = sqrt(dot_self(z[j]-z[i])) + sqrt(dot_self(z[j]-z_s[i]));
    }
  }
}

model {
  for (i in 1:M){
    z_s[i] ~ normal(0,rho_s);
    z[i] ~ normal(0, rho);
  }
  {
    vector[K] f;
    real eta_ij;
    for (i in 1:M) {
      for (j in 1:N) {
        if (!w[i, j]) {
          eta_ij = D[i,j] + asym[i, j];
          f[1] = 1 - Phi( (eta_ij - b[1]) / sigma);
          for(k in 2:(K-1))
            f[k] = Phi( (eta_ij - b[k-1]) / sigma) - Phi( (eta_ij - b[k]) / sigma);
          f[K] = Phi( (eta_ij - b[K-1] ) / sigma);
          mat[i, j] ~ categorical(f);
        }
      }
    }
  }
}
