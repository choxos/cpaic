// Component-additive ML-NMR for exact survival outcomes with a continuous
// M-spline baseline. The hazard basis, integrated basis, censoring and
// delayed-entry likelihoods, and likelihood-level AgD integration are ported
// from multinma (Phillippo et al. 2020). Both packages are GPL-3.

functions {
  real cpaic_survival_loglik(
      row_vector time_basis,
      row_vector itime_basis,
      row_vector start_basis,
      row_vector entry_basis,
      int delayed,
      int status,
      real eta,
      vector coefficients) {
    // Every contribution is conditioned on survival to the entry time a (left
    // truncation), so each is written through the ENTRY-CONDITIONED cumulative
    // hazard H(t) - H(a). The conditioning must be applied INSIDE the censoring
    // probability, not added afterwards: (1 - S(t)) / S(a) is not a probability
    // and can exceed one, whereas 1 - S(t)/S(a) is the correct left-censoring
    // term.
    real entry_cum = dot_product(entry_basis, coefficients) * exp(eta);
    real cum = dot_product(itime_basis, coefficients) * exp(eta) - entry_cum;
    real out;

    if (status == 0) {
      out = -cum;                                     // S(t) / S(a)
    } else if (status == 1) {
      out = -cum + log(dot_product(time_basis, coefficients)) + eta;
    } else if (status == 2) {
      out = log1m_exp(-cum);                          // 1 - S(t) / S(a)
    } else {
      real start_cum =
        dot_product(start_basis, coefficients) * exp(eta) - entry_cum;
      out = log_diff_exp(-start_cum, -cum);           // (S(s) - S(t)) / S(a)
    }
    return out;
  }

  real cpaic_event_probability(
      row_vector itime_basis,
      row_vector entry_basis,
      real eta,
      vector coefficients) {
    real cumulative = (dot_product(itime_basis, coefficients)
                       - dot_product(entry_basis, coefficients)) * exp(eta);
    return -expm1(-fmax(cumulative, 0));
  }
}

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=1> P;
  int<lower=1> Q;
  int<lower=1> n_int;
  int<lower=4> N_base;

  matrix[N_ipd, N_base] time_basis_ipd;
  matrix[N_ipd, N_base] itime_basis_ipd;
  matrix[N_ipd, N_base] start_basis_ipd;
  matrix[N_ipd, N_base] entry_basis_ipd;
  array[N_ipd] int<lower=0, upper=1> delayed_ipd;
  array[N_ipd] int<lower=0, upper=3> status_ipd;
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;
  matrix[N_ipd, P] X_ipd;
  array[N_ipd, Q] int em_idx;

  matrix[N_agd, N_base] time_basis_agd;
  matrix[N_agd, N_base] itime_basis_agd;
  matrix[N_agd, N_base] start_basis_agd;
  matrix[N_agd, N_base] entry_basis_agd;
  array[N_agd] int<lower=0, upper=1> delayed_agd;
  array[N_agd] int<lower=0, upper=3> status_agd;
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
  simplex[N_base] coefficients;
  vector[C] beta;
  vector[P] breg;
  matrix[C, Q] gamma;
  vector[RE ? N_delta : 0] delta_aux;
  vector<lower=0>[RE] tau;
}

transformed parameters {
  vector[RE ? N_delta : 0] delta;
  vector[N_ipd] eta_ipd;
  vector[N_ipd] log_L_ipd;
  vector[N_agd] log_L_agd;
  vector[N_ipd] p_event_ipd;
  vector[N_agd] p_event_agd;

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
    log_L_ipd[i] = cpaic_survival_loglik(
      time_basis_ipd[i], itime_basis_ipd[i], start_basis_ipd[i],
      entry_basis_ipd[i], delayed_ipd[i], status_ipd[i], lp, coefficients
    );
    p_event_ipd[i] = cpaic_event_probability(
      itime_basis_ipd[i], entry_basis_ipd[i], lp, coefficients
    );
  }

  for (a in 1:N_agd) {
    vector[n_int] log_L_ii;
    real event_acc = 0;
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
      log_L_ii[k] = cpaic_survival_loglik(
        time_basis_agd[a], itime_basis_agd[a], start_basis_agd[a],
        entry_basis_agd[a], delayed_agd[a], status_agd[a], lp, coefficients
      );
      event_acc += cpaic_event_probability(
        itime_basis_agd[a], entry_basis_agd[a], lp, coefficients
      );
    }
    log_L_agd[a] = log_sum_exp(log_L_ii) - log(n_int);
    p_event_agd[a] = event_acc / n_int;
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  coefficients ~ dirichlet(rep_vector(1.0, N_base));
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
    target += sum(log_L_ipd);
    target += sum(log_L_agd);
  }
}

generated quantities {
  vector[N_ipd + N_agd] log_lik;
  array[N_ipd] int event_rep_ipd;
  array[N_agd] int event_rep_agd;

  for (i in 1:N_ipd) {
    log_lik[i] = log_L_ipd[i];
    event_rep_ipd[i] = bernoulli_rng(fmin(fmax(p_event_ipd[i], 0), 1));
  }
  for (a in 1:N_agd) {
    log_lik[N_ipd + a] = log_L_agd[a];
    event_rep_agd[a] = bernoulli_rng(fmin(fmax(p_event_agd[a], 0), 1));
  }
}
