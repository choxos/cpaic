# Construct a disconnected component network with IPD and known truth.
#
# Structure (binary outcome, log-OR scale):
#   Sub-network 1 (anchored on Placebo): Placebo, A, B
#     S1 (AgD): A      vs Placebo   -> beta_A
#     S2 (AgD): B      vs Placebo   -> beta_B
#   Sub-network 2 (isolated; bridged via shared components A, B):
#     S3 (IPD): A+B+C  vs A+B       -> beta_C (+ effect modification)
#     S4 (IPD): A+B+D  vs A+B       -> beta_D (+ effect modification)
#     S5 (AgD): A+B+C  vs A+B+D     -> beta_C - beta_D
#
# No treatment is shared between the two sub-networks, so the network is
# disconnected; components A and B bridge it, and C, D are identified
# within sub-network 2. One effect modifier x1 is imbalanced in the IPD
# studies, so population adjustment matters.

set.seed(2026)

# True component effects (log-OR) and effect modification on x1.
beta <- c(A = 0.5, B = 0.4, C = 0.3, D = 0.6)
em   <- c(C = 0.8, D = -0.7)   # treatment x x1 interaction for C and D arms
prog_x1 <- 0.5                  # prognostic effect of x1
intercept <- -0.5

# ---- IPD studies S3 and S4 (x1 shifted away from the target mean 0) ----
gen_ipd <- function(study, active_arm, ref_arm, comp, n = 1600, mu_x1 = 1.4) {
  x1 <- stats::rnorm(n, mean = mu_x1, sd = 1)
  arm <- rep(c(active_arm, ref_arm), each = n / 2)
  is_active <- arm == active_arm
  eta <- intercept + beta[[comp]] * is_active + prog_x1 * x1 +
    em[[comp]] * is_active * x1
  y <- stats::rbinom(n, 1, plogis(eta))
  data.frame(.study = study, .trt = arm, .y = y, x1 = x1,
             stringsAsFactors = FALSE)
}

ipd_S3 <- gen_ipd("S3", "A+B+C", "A+B", "C", mu_x1 = 1.5)
ipd_S4 <- gen_ipd("S4", "A+B+D", "A+B", "D", mu_x1 = 1.3)
cpaic_bin_ipd <- rbind(ipd_S3, ipd_S4)

# Naive (unadjusted) contrasts for the IPD studies, so the AgD network is
# complete even before any population adjustment.
naive_contrast <- function(d, active, ref) {
  d$.arm <- stats::relevel(factor(d$.trt), ref = ref)
  g <- stats::glm(.y ~ .arm, family = stats::binomial(), data = d)
  s <- summary(g)$coefficients[2, 1:2]
  c(TE = unname(s[1]), seTE = unname(s[2]))
}
c3 <- naive_contrast(ipd_S3, "A+B+C", "A+B")
c4 <- naive_contrast(ipd_S4, "A+B+D", "A+B")

# ---- Aggregate contrasts ----
cpaic_bin_agd <- data.frame(
  studlab = c("S1", "S2", "S3", "S4", "S5"),
  treat1  = c("A", "B", "A+B+C", "A+B+D", "A+B+C"),
  treat2  = c("Placebo", "Placebo", "A+B", "A+B", "A+B+D"),
  TE      = c(beta[["A"]], beta[["B"]], c3[["TE"]], c4[["TE"]],
              beta[["C"]] - beta[["D"]]),
  seTE    = c(0.12, 0.12, c3[["seTE"]], c4[["seTE"]], 0.15),
  stringsAsFactors = FALSE
)

attr(cpaic_bin_agd, "truth") <- beta
attr(cpaic_bin_ipd, "truth") <- list(beta = beta, em = em,
                                     prog_x1 = prog_x1,
                                     target_x1 = 0)

usethis::use_data(cpaic_bin_agd, cpaic_bin_ipd, overwrite = TRUE)
