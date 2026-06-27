// Component-additive ML-NMR, continuous (normal) outcome, identity link.
// With an identity link the population-average mean is linear in the
// covariates, so a single integration point at the covariate means is
// exact; the model nevertheless averages over the supplied points so the
// same data interface is used for every family.

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=1> n_int;

  vector[N_ipd] y_ipd;
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;
  matrix[N_ipd, P] X_ipd;
  array[N_ipd, Q] int em_idx;

  vector[N_agd] y_agd;           // arm mean outcome
  vector<lower=0>[N_agd] se_agd; // known standard error of the mean
  array[N_agd] int<lower=1, upper=N_studies> study_agd;
  matrix[N_agd, C] Tc_agd;
  matrix[N_agd * n_int, P] X_agd_int;

  real prior_intercept_sd;
  real prior_beta_sd;
  real prior_reg_sd;
}

transformed data {
  array[Q] int emc;
  for (q in 1:Q) emc[q] = em_idx[1, q];
}

parameters {
  vector[N_studies] mu;
  vector[C] beta;
  vector[P] breg;
  matrix[C, Q] gamma;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N_ipd] eta_ipd;
  vector[N_agd] eta_agd;

  for (i in 1:N_ipd) {
    real lp = mu[study_ipd[i]] + X_ipd[i] * breg + Tc_ipd[i] * beta;
    for (q in 1:Q)
      lp += dot_product(Tc_ipd[i], gamma[, q]) * X_ipd[i, emc[q]];
    eta_ipd[i] = lp;
  }
  for (a in 1:N_agd) {
    real acc = 0;
    for (k in 1:n_int) {
      int row = (a - 1) * n_int + k;
      real lp = mu[study_agd[a]] + X_agd_int[row] * breg + Tc_agd[a] * beta;
      for (q in 1:Q)
        lp += dot_product(Tc_agd[a], gamma[, q]) * X_agd_int[row, emc[q]];
      acc += lp;
    }
    eta_agd[a] = acc / n_int;
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  to_vector(gamma) ~ normal(0, prior_reg_sd);
  sigma ~ normal(0, prior_beta_sd);

  if (N_ipd > 0) y_ipd ~ normal(eta_ipd, sigma);
  if (N_agd > 0) y_agd ~ normal(eta_agd, se_agd);
}
