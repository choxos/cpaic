// Component-additive ML-NMR, survival outcome, piecewise-exponential (PWE)
// baseline hazard. The baseline log-hazard is a step function with K
// interval-specific levels (bshape), so K = 1 recovers the exponential
// model and larger K gives an arbitrarily flexible baseline. The hazard is
// h(t | x) = exp(mu_study + bshape[interval(t)] + x'breg + Tc'beta + ...),
// proportional across treatments (component-additive log hazard ratios).
//
// IPD are supplied Lexis-expanded: one row per individual-interval at risk,
// with the time at risk in that interval and a 0/1 event indicator. The
// exponential likelihood contribution then equals a Poisson log-likelihood
// for the event indicator with log(time-at-risk) offset. Aggregate arms are
// supplied per interval (events and person-time) and fitted by integrating
// the hazard over the covariate distribution.

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=1> n_int;
  int<lower=1> K;                              // number of baseline intervals

  array[N_ipd] int<lower=0, upper=1> y_ipd;    // event in this interval-row
  vector<lower=0>[N_ipd] time_ipd;             // time at risk in interval
  array[N_ipd] int<lower=1, upper=K> interval_ipd;
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;
  matrix[N_ipd, P] X_ipd;
  array[N_ipd, Q] int em_idx;

  array[N_agd] int<lower=0> r_agd;             // events in arm-interval
  vector<lower=0>[N_agd] E_agd;                // person-time in arm-interval
  array[N_agd] int<lower=1, upper=K> interval_agd;
  array[N_agd] int<lower=1, upper=N_studies> study_agd;
  matrix[N_agd, C] Tc_agd;
  matrix[N_agd * n_int, P] X_agd_int;

  real prior_intercept_sd;
  real prior_beta_sd;
  real prior_reg_sd;
}

transformed data {
  array[Q] int emc;
  vector[N_ipd] log_time;
  for (q in 1:Q) emc[q] = em_idx[1, q];
  for (i in 1:N_ipd) log_time[i] = log(time_ipd[i]);
}

parameters {
  vector[N_studies] mu;
  vector[K - 1] bshape_raw;       // baseline shape; interval 1 fixed at 0
  vector[C] beta;
  vector[P] breg;
  matrix[C, Q] gamma;
}

transformed parameters {
  vector[K] bshape = append_row(rep_vector(0, 1), bshape_raw);
  vector[N_ipd] eta_ipd;
  vector[N_agd] lambda_agd;

  for (i in 1:N_ipd) {
    real lp = mu[study_ipd[i]] + bshape[interval_ipd[i]] + X_ipd[i] * breg
              + Tc_ipd[i] * beta;
    for (q in 1:Q)
      lp += dot_product(Tc_ipd[i], gamma[, q]) * X_ipd[i, emc[q]];
    eta_ipd[i] = lp + log_time[i];
  }
  for (a in 1:N_agd) {
    real acc = 0;
    for (k in 1:n_int) {
      int row = (a - 1) * n_int + k;
      real lp = mu[study_agd[a]] + bshape[interval_agd[a]]
                + X_agd_int[row] * breg + Tc_agd[a] * beta;
      for (q in 1:Q)
        lp += dot_product(Tc_agd[a], gamma[, q]) * X_agd_int[row, emc[q]];
      acc += exp(lp);
    }
    lambda_agd[a] = E_agd[a] * acc / n_int;
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  bshape_raw ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  to_vector(gamma) ~ normal(0, prior_reg_sd);

  if (N_ipd > 0) y_ipd ~ poisson_log(eta_ipd);
  if (N_agd > 0) r_agd ~ poisson(lambda_agd);
}
