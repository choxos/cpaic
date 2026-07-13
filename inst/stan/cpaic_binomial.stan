// Component-additive ML-NMR for binary outcomes with a logit link.
// Random-effects and integration logic follow multinma (Phillippo et al.
// 2020), reparameterized for component effects. Both packages are GPL-3.

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=1> P;
  int<lower=1> Q;
  int<lower=1> n_int;

  array[N_ipd] int<lower=0, upper=1> y_ipd;
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;
  matrix[N_ipd, P] X_ipd;
  array[N_ipd, Q] int em_idx;

  array[N_agd] int<lower=0> r_agd;
  array[N_agd] int<lower=1> n_agd;
  array[N_agd] int<lower=1, upper=N_studies> study_agd;
  matrix[N_agd, C] Tc_agd;
  matrix[N_agd * n_int, P] X_agd_int;

  int<lower=0, upper=1> RE;
  int<lower=0, upper=1> noncentered;
  int<lower=1> N_delta;
  matrix[N_delta, N_delta] L_delta;
  array[N_ipd] int<lower=0, upper=N_delta> re_idx_ipd;
  array[N_agd] int<lower=0, upper=N_delta> re_idx_agd;

  int<lower=0, upper=1> prior_only;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_reg_sd;
  int<lower=1, upper=2> prior_gamma_dist;
  real<lower=0> prior_gamma_scale;
  real<lower=0> prior_gamma_df;
  int<lower=1, upper=2> prior_tau_dist;
  real<lower=0> prior_tau_scale;
  real<lower=0> prior_tau_df;
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
  vector[RE ? N_delta : 0] delta_aux;
  vector<lower=0>[RE] tau;
}

transformed parameters {
  vector[RE ? N_delta : 0] delta;
  vector[N_ipd] eta_ipd;
  vector[N_agd] p_agd;

  delta = delta_aux;
  if (RE) {
    if (noncentered) delta = tau[1] * L_delta * delta_aux;
  }

  for (i in 1:N_ipd) {
    real lp = mu[study_ipd[i]] + X_ipd[i] * breg + Tc_ipd[i] * beta;
    for (q in 1:Q) {
      lp += dot_product(Tc_ipd[i], gamma[, q]) * X_ipd[i, emc[q]];
    }
    if (RE) {
      if (re_idx_ipd[i] > 0) lp += delta[re_idx_ipd[i]];
    }
    eta_ipd[i] = lp;
  }

  for (a in 1:N_agd) {
    real acc = 0;
    for (k in 1:n_int) {
      int row = (a - 1) * n_int + k;
      real lp = mu[study_agd[a]] + X_agd_int[row] * breg
                + Tc_agd[a] * beta;
      for (q in 1:Q) {
        lp += dot_product(Tc_agd[a], gamma[, q]) * X_agd_int[row, emc[q]];
      }
      if (RE) {
        if (re_idx_agd[a] > 0) lp += delta[re_idx_agd[a]];
      }
      acc += inv_logit(lp);
    }
    p_agd[a] = acc / n_int;
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  if (prior_gamma_dist == 1) {
    to_vector(gamma) ~ normal(0, prior_gamma_scale);
  } else {
    to_vector(gamma) ~ student_t(prior_gamma_df, 0, prior_gamma_scale);
  }

  if (RE) {
    if (prior_tau_dist == 1) tau ~ normal(0, prior_tau_scale);
    else tau ~ student_t(prior_tau_df, 0, prior_tau_scale);
    if (noncentered) {
      delta_aux ~ std_normal();
    } else {
      delta_aux ~ multi_normal_cholesky(
        rep_vector(0, N_delta), tau[1] * L_delta
      );
    }
  }

  if (!prior_only) {
    if (N_ipd > 0) y_ipd ~ bernoulli_logit(eta_ipd);
    if (N_agd > 0) r_agd ~ binomial(n_agd, p_agd);
  }
}

generated quantities {
  vector[N_ipd + N_agd] log_lik;
  array[N_ipd] int yrep_ipd;
  array[N_agd] int rrep_agd;

  for (i in 1:N_ipd) {
    log_lik[i] = bernoulli_logit_lpmf(y_ipd[i] | eta_ipd[i]);
    yrep_ipd[i] = bernoulli_logit_rng(eta_ipd[i]);
  }
  for (a in 1:N_agd) {
    log_lik[N_ipd + a] = binomial_lpmf(r_agd[a] | n_agd[a], p_agd[a]);
    rrep_agd[a] = binomial_rng(n_agd[a], p_agd[a]);
  }
}
