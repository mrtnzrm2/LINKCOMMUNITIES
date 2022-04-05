functions {
  int[] find_interval(vector w, real inter_1, real inter_2, int N) {
    int interval[2] = {0, 0};
    for (i in 1:N) {
      if (w[i] <= inter_1){
        interval[1] = i;
      }
      if (w[i] <= inter_2)
        interval[2] = i;
    }
    return interval;
  }
  real dispersion_interval(vector dist, int[] intervals, int N) {
    vector[intervals[2] - intervals[1] + 1] split_dist;
    int count = 1;
    for (i in 1:N) {
      if (i >= intervals[1] && i <= intervals[2]) {
        split_dist[count] = dist[i];
        count = count + 1;
      }
    }
    return sd(split_dist);
  }
}

data {
  int<lower=2, upper=15> K;
  int<lower=100> n;
  vector[2] range;
  vector[n] w;
  vector[n] dist;
}

parameters {
  ordered[K - 1] intervals;
}

model {
  intervals ~ uniform(range[1], range[2]);
  {
    int inter[2];
    vector[K] disp_interval;
    for (i in 1:K) {
      if (i == 1)
        inter = find_interval(
          w, range[1], intervals[1], n
        );
      else if (i > 1 && i < K - 1)
        inter = find_interval(
          w, intervals[i], intervals[i + 1], n
        );
      else
        inter = find_interval(
          w, intervals[K - 1], range[2], n
        );
      disp_interval[i] = dispersion_interval(dist, inter, n);
    }
    disp_interval ~ normal(0, 0.01);
  }
}
