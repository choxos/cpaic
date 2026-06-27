# Package index

## Network setup

Build a (possibly disconnected) component network and code components.

- [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  : Set up a (possibly disconnected) component network for cpaic
- [`build_C_matrix()`](https://choxos.github.io/cpaic/reference/build_C_matrix.md)
  : Build a component-coded treatment-by-component matrix

## Connectivity

Detect disconnection and check component-bridge identifiability.

- [`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
  : Assess connectivity and component-bridge identifiability of a
  network

## Connection layer

Reconnect a disconnected network through shared components.

- [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
  : Reconnect a network through its additive component structure
- [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md)
  : Test the additivity assumption of the component model

## Population adjustment

Anchored, population-adjusted indirect comparison across the network.

- [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) :
  Component matching-adjusted indirect comparison (cMAIC)
- [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) :
  Component simulated treatment comparison (cSTC)
- [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) :
  Component-additive multilevel network meta-regression (ML-NMR)
- [`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md)
  : Effective sample sizes from a cMAIC fit

## Reporting and visualization

- [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
  : Relative treatment effects from a cpaic fit
- [`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md)
  : League table of all pairwise relative effects
- [`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md)
  : Component effects from a cpaic fit
- [`forest()`](https://choxos.github.io/cpaic/reference/forest.md) :
  Forest plot of relative effects
- [`plot(`*`<cpaic_network>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_network.md)
  : Plot the component network

## Data

- [`cpaic_bin_agd`](https://choxos.github.io/cpaic/reference/cpaic_bin_agd.md)
  : Example disconnected component network: aggregate contrasts
- [`cpaic_bin_ipd`](https://choxos.github.io/cpaic/reference/cpaic_bin_ipd.md)
  : Example disconnected component network: individual patient data
