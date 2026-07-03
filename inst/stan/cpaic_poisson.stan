// Component-additive ML-NMR, count outcome, Poisson / log link. Also used
// for exponential survival (events with person-time exposure as offset).
// Aggregate arms are fitted by integrating the individual rate over the
// covariate distribution: expected count = exposure * mean(exp(eta)).

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=1> n_int;

  array[N_ipd] int<lower=0> y_ipd;        // counts (or 0/1 event indicator)
  vector<lower=0>[N_ipd] offset_ipd;      // log already applied outside? no:
                                          // raw exposure; log taken here
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;
  matrix[N_ipd, P] X_ipd;
  array[N_ipd, Q] int em_idx;

  array[N_agd] int<lower=0> r_agd;        // total events
  vector<lower=0>[N_agd] E_agd;           // total exposure (person-time)
  array[N_agd] int<lower=1, upper=N_studies> study_agd;
  matrix[N_agd, C] Tc_agd;
  matrix[N_agd * n_int, P] X_agd_int;

  real prior_intercept_sd;
  real prior_beta_sd;
  real prior_reg_sd;
}

transformed data {
  array[Q] int emc;
  vector[N_ipd] log_offset;
  for (q in 1:Q) emc[q] = em_idx[1, q];
  for (i in 1:N_ipd) log_offset[i] = log(offset_ipd[i]);
}

parameters {
  vector[N_studies] mu;
  vector[C] beta;
  vector[P] breg;
  matrix[C, Q] gamma;
}

transformed parameters {
  vector[N_ipd] eta_ipd;
  vector[N_agd] lambda_agd;       // expected total events per aggregate arm

  for (i in 1:N_ipd) {
    real lp = mu[study_ipd[i]] + X_ipd[i] * breg + Tc_ipd[i] * beta;
    for (q in 1:Q)
      lp += dot_product(Tc_ipd[i], gamma[, q]) * X_ipd[i, emc[q]];
    eta_ipd[i] = lp + log_offset[i];
  }
  for (a in 1:N_agd) {
    real acc = 0;
    for (k in 1:n_int) {
      int row = (a - 1) * n_int + k;
      real lp = mu[study_agd[a]] + X_agd_int[row] * breg + Tc_agd[a] * beta;
      for (q in 1:Q)
        lp += dot_product(Tc_agd[a], gamma[, q]) * X_agd_int[row, emc[q]];
      acc += exp(lp);
    }
    lambda_agd[a] = E_agd[a] * acc / n_int;   // exposure * mean rate
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  to_vector(gamma) ~ normal(0, prior_reg_sd);

  if (N_ipd > 0) y_ipd ~ poisson_log(eta_ipd);
  if (N_agd > 0) r_agd ~ poisson(lambda_agd);
}
