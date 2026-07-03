# Plot the component network

Draws the treatment network, coloring nodes by sub-network so that
disconnection is visible at a glance. Edges with population-adjusted IPD
evidence are not distinguished here (see the vignettes).

## Usage

``` r
# S3 method for class 'cpaic_network'
plot(x, ...)
```

## Arguments

- x:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  object.

- ...:

  Passed to
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

## Value

The `igraph` object, invisibly.
