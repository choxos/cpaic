# cpaic: Component-Based Population-Adjusted Indirect Comparison

Experimental research software that combines component network
meta-analysis (cNMA) with population-adjusted indirect comparison (PAIC)
methods for disconnected treatment networks. The additive component
structure of cNMA (Ruecker et al. 2020
[doi:10.1002/bimj.201800167](https://doi.org/10.1002/bimj.201800167) )
can bridge sub-networks that share treatment components. Anchored
simulated treatment comparison (STC) and matching-adjusted indirect
comparison (MAIC) adjust individual-patient-data contrasts before a
two-stage component bridge. This can mix target-specific contrasts with
retained aggregate contrasts from different study populations; marginal
contrasts on nonlinear scales also need not be component-additive, so
affected bridges require explicit experimental opt-in.
Component-additive multilevel network meta-regression (ML-NMR; Phillippo
et al. 2020 [doi:10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579)
) jointly models individual and aggregate data. Its target summaries are
average conditional effects on the link scale evaluated at supplied
covariate means, not marginally standardized odds, rate, or hazard
ratios. Models are implemented for binary, continuous, count, and
time-to-event outcomes. This methodology and implementation have not
been validated for clinical, regulatory, reimbursement, or other
decision use.

## See also

Useful links:

- <https://github.com/choxos/cpaic>

- <https://choxos.github.io/cpaic/>

- Report bugs at <https://github.com/choxos/cpaic/issues>

## Author

**Maintainer**: Ahmad Sofi-Mahmudi <ahmad.pub@gmail.com>
([ORCID](https://orcid.org/0000-0001-6829-0823))

Authors:

- Ahmad Sofi-Mahmudi <ahmad.pub@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-6829-0823))
