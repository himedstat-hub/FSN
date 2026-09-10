data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> G;
  matrix[N, K] X;
  array[N] int<lower=1, upper=G> diagnosis;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma_dx;
  vector[G] z_dx;
}
transformed parameters {
  vector[G] dx_effect = z_dx * sigma_dx;
}
model {
  alpha ~ normal(0, 2.5);
  beta ~ normal(0, 1);
  sigma_dx ~ normal(0, 0.5);
  z_dx ~ normal(0, 1);
  y ~ bernoulli_logit(alpha + X * beta + dx_effect[diagnosis]);
}
generated quantities {
  vector[N] p;
  for (n in 1:N) {
    p[n] = inv_logit(alpha + X[n] * beta + dx_effect[diagnosis[n]]);
  }
}

