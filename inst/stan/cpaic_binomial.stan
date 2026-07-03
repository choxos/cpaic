// Component-additive multilevel network meta-regression (ML-NMR), binary
// outcome, logit link. The treatment effect of an arm is the sum of its
// component effects (design row Tc times beta), so disconnected
// sub-networks that share components are connected through shared
// component parameters. Aggregate-data arms are fitted by numerically
// integrating the individual-level model over the covariate distribution
// of the aggregate population (integration points supplied as data).
//
// Adapted in spirit from the ML-NMR implementation in the 'multinma'
// package (Phillippo et al. 2020), reparameterized to component effects.

data {
  // Dimensions
  int<lower=0> N_ipd;            // individuals with IPD
  int<lower=0> N_agd;            // aggregate arms
  int<lower=1> N_studies;        // total studies (IPD + AgD)
  int<lower=1> C;                // number of components
  int<lower=0> P;                // number of covariates (prognostic + EM)
  int<lower=0> Q;                // number of effect modifiers (subset of P)
  int<lower=1> n_int;            // integration points per AgD arm

  // IPD block
  array[N_ipd] int<lower=0, upper=1> y_ipd;
  array[N_ipd] int<lower=1, upper=N_studies> study_ipd;
  matrix[N_ipd, C] Tc_ipd;       // component design (arm membership)
  matrix[N_ipd, P] X_ipd;        // covariates
  array[N_ipd, Q] int em_idx;    // columns of X that are effect modifiers
                                 // (1-based; same for all rows)

  // AgD block (integration points stacked: row (a-1)*n_int + q)
  array[N_agd] int<lower=0> r_agd;
  array[N_agd] int<lower=1> n_agd;
  array[N_agd] int<lower=1, upper=N_studies> study_agd;
  matrix[N_agd, C] Tc_agd;
  matrix[N_agd * n_int, P] X_agd_int;

  // Priors
  real prior_intercept_sd;
  real prior_beta_sd;
  real prior_reg_sd;
}

transformed data {
  array[Q] int emc;
  for (q in 1:Q) emc[q] = em_idx[1, q];
}

parameters {
  vector[N_studies] mu;          // study intercepts (baselines)
  vector[C] beta;                // component main effects
  vector[P] breg;                // covariate main effects (prognostic)
  matrix[C, Q] gamma;            // component x effect-modifier interactions
}

transformed parameters {
  vector[N_ipd] eta_ipd;
  vector[N_agd] p_agd;

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
      acc += inv_logit(lp);
    }
    p_agd[a] = acc / n_int;       // population-average event probability
  }
}

model {
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
  breg ~ normal(0, prior_reg_sd);
  to_vector(gamma) ~ normal(0, prior_reg_sd);

  if (N_ipd > 0) y_ipd ~ bernoulli_logit(eta_ipd);
  if (N_agd > 0) r_agd ~ binomial(n_agd, p_agd);
}

generated quantities {
  // Component effects are `beta`; treatment effects are recovered as C %*% beta
  // on the R side from the component design matrix.
}
