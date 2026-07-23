// Component-additive ML-NMR for continuous outcomes with an identity link.
// Random-effects, integration logic, and the scaled thin QR approach follow
// multinma (Phillippo et al. 2020), reparameterized for component effects.
// Both packages are GPL-3.

data {
  int<lower=0> N_ipd;
  int<lower=0> N_agd;
  int<lower=1> N_studies;
  int<lower=1> C;
  int<lower=1> P;
  int<lower=1> Q;
  int<lower=1> n_int;
  int<lower=1> nX;

  vector[N_ipd] y_ipd;
  matrix[N_ipd, nX] Z_ipd;

  vector[N_agd] y_agd;
  vector<lower=0>[N_agd] se_agd;
  matrix[N_agd * n_int, nX] Z_agd_int;

  int<lower=0, upper=1> QR;
  matrix[QR ? nX : 0, QR ? nX : 0] R_inv;

  int<lower=0, upper=1> RE;
  int<lower=0, upper=1> noncentered;
  int<lower=1> N_delta;
  matrix[N_delta, N_delta] L_delta;
  array[N_ipd] int<lower=0, upper=N_delta> re_idx_ipd;
  array[N_agd] int<lower=0, upper=N_delta> re_idx_agd;

  int<lower=0, upper=1> prior_only;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_reg_sd;
  int<lower=1, upper=2> prior_gamma_dist;
  real<lower=0> prior_gamma_scale;
  real<lower=0> prior_gamma_df;
  int<lower=1, upper=2> prior_tau_dist;
  real<lower=0> prior_tau_scale;
  real<lower=0> prior_tau_df;
}

parameters {
  vector[nX] beta_tilde;
  real<lower=0> sigma;
  vector[RE ? N_delta : 0] delta_aux;
  vector<lower=0>[RE] tau;
}

transformed parameters {
  vector[nX] allbeta;
  vector[N_studies] mu;
  vector[C] beta;
  vector[P] breg;
  matrix[C, Q] gamma;
  vector[RE ? N_delta : 0] delta;
  vector[N_ipd] eta_ipd;
  vector[N_agd] eta_agd;

  allbeta = QR ? R_inv * beta_tilde : beta_tilde;
  mu = segment(allbeta, 1, N_studies);
  beta = segment(allbeta, N_studies + 1, C);
  breg = segment(allbeta, N_studies + C + 1, P);
  gamma = to_matrix(
    segment(allbeta, N_studies + C + P + 1, C * Q), C, Q
  );

  delta = delta_aux;
  if (RE) {
    if (noncentered) delta = tau[1] * L_delta * delta_aux;
  }

  for (i in 1:N_ipd) {
    real lp = QR ? Z_ipd[i] * beta_tilde : Z_ipd[i] * allbeta;
    if (RE) {
      if (re_idx_ipd[i] > 0) lp += delta[re_idx_ipd[i]];
    }
    eta_ipd[i] = lp;
  }
  for (a in 1:N_agd) {
    real acc = 0;
    for (k in 1:n_int) {
      int row = (a - 1) * n_int + k;
      real lp = QR ? Z_agd_int[row] * beta_tilde
                   : Z_agd_int[row] * allbeta;
      if (RE) {
        if (re_idx_agd[a] > 0) lp += delta[re_idx_agd[a]];
      }
      acc += lp;
    }
    eta_agd[a] = acc / n_int;
  }
}

model {
  // Priors remain on allbeta. The QR map is linear with a constant Jacobian,
  // so no Jacobian adjustment can change the posterior distribution.
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  if (prior_gamma_dist == 1) {
    to_vector(gamma) ~ normal(0, prior_gamma_scale);
  } else {
    to_vector(gamma) ~ student_t(prior_gamma_df, 0, prior_gamma_scale);
  }
  sigma ~ normal(0, prior_sigma_sd);   // residual SD: its own prior scale

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
    if (N_ipd > 0) y_ipd ~ normal(eta_ipd, sigma);
    if (N_agd > 0) y_agd ~ normal(eta_agd, se_agd);
  }
}

generated quantities {
  vector[N_ipd + N_agd] log_lik;
  vector[N_ipd] yrep_ipd;
  vector[N_agd] yrep_agd;

  for (i in 1:N_ipd) {
    log_lik[i] = normal_lpdf(y_ipd[i] | eta_ipd[i], sigma);
    yrep_ipd[i] = normal_rng(eta_ipd[i], sigma);
  }
  for (a in 1:N_agd) {
    log_lik[N_ipd + a] = normal_lpdf(y_agd[a] | eta_agd[a], se_agd[a]);
    yrep_agd[a] = normal_rng(eta_agd[a], se_agd[a]);
  }
}
