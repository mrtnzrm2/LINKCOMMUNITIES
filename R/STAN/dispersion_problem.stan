functions {
  array[] int find_interval(vector w, real inter_1, real inter_2, int N) {
    array[2] int interval = {1, 1};
    for (i in 1:N) {
      if (w[i] <= inter_1){
        interval[1] = i;
      }
      if (w[i] <= inter_2)
        interval[2] = i;
    }
    return interval;
  }
  real dispersion_interval(vector dist, array[] int intervals, int N) {
    int size_vector  = intervals[2] - intervals[1] + 1; 
    vector[size_vector] split_dist;
    real sd_dist;
    int count = 1;
    for (i in 1:N) {
      if (i >= intervals[1] && i <= intervals[2]) {
        split_dist[count] = dist[i];
        count = count + 1;
      }
    }
    if (is_nan(split_dist[size_vector])) {
      if (size_vector > 2)
        sd_dist = sd(split_dist[1:(size_vector - 1)]);
      else 
        sd_dist = 0;
    }
    else
      sd_dist = sd(split_dist);
    return sd_dist;
  }
}

data {
  int<lower=2, upper=15> K;
  int<lower=100> n;
  vector[2] range;
  vector[n] w;
  vector[n] dist;
}

transformed data {
  real delta = 1e-3 + 1e-5;
}

parameters {
  ordered[K - 1] intervals;
}

transformed parameters {
  vector[K + 1] Intervals;
  Intervals[1] = range[1];
  Intervals[K + 1] = range[2];
  for (i in 1:(K - 1))
    Intervals[i + 1] = intervals[i];
}

model {
  intervals ~ uniform(range[1] + delta, range[2] - delta);
  {
    array [2] int inter;
    vector[K] disp_interval;
    for (i in 1:K) {
      inter = find_interval(
          w, Intervals[i], Intervals[i + 1], n
      );
      disp_interval[i] = dispersion_interval(dist, inter, n);
    }
    disp_interval ~ normal(0, 0.01);
  }
}
