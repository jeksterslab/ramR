Reticular Action Model (RAM) Notation
================
Ivan Jacob Agaloos Pesigan
2021-01-23

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R build
status](https://github.com/jeksterslab/ramR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/jeksterslab/ramR/actions?workflow=R-CMD-check)
[![Travis build
status](https://travis-ci.com/jeksterslab/ramR.svg?branch=master)](https://travis-ci.com/jeksterslab/ramR)
[![codecov](https://codecov.io/github/jeksterslab/ramR/branch/master/graphs/badge.svg)](https://codecov.io/github/jeksterslab/ramR)
<!-- badges: end -->

## Description

A collection of utility functions using the Reticular Action Model (RAM)
notation.

## Installation

You can install the released version of `ramR` from
[GitHub](https://github.com/jeksterslab/ramR) with:

``` r
remotes::install_github("jeksterslab/ramR")
```

## Symbolic Example

This is a symbolic example for the model

*y* = *α* + *β* ⋅ *x* + *ε*.

``` r
A <- S <- matrix(
  data = 0,
  nrow = 3,
  ncol = 3
)
A[1, ] <- c(0, "beta", 1)
diag(S) <- c(0, "sigma[x]^2", "sigma[varepsilon]^2")
filter <- diag(2)
filter <- cbind(filter, 0)
u <- c("alpha", "mu[x]", 0)
```

The covariance expectations can be symbolically derived using the
`ramR::C_sym()` function.

``` r
ramR::C_sym(A, S)
#> {{sigma[x]^2*beta^2+sigma[varepsilon]^2,                       beta*sigma[x]^2,                   sigma[varepsilon]^2},
#>  {                      sigma[x]^2*beta,                            sigma[x]^2,                                     0},
#>  {                  sigma[varepsilon]^2,                                     0,                   sigma[varepsilon]^2}}
```

The covariance expectations for the observed variables can be
symbolically derived using the `ramR::M_sym()` function.

``` r
ramR::M_sym(A, S, filter)
#> {{sigma[x]^2*beta^2+sigma[varepsilon]^2,                       beta*sigma[x]^2},
#>  {                      sigma[x]^2*beta,                            sigma[x]^2}}
```

The mean expectations can be symbolically derived using the
`ramR::v_sym()` function.

``` r
ramR::v_sym(A, u)
#> {{alpha+beta*mu[x]},
#>  {           mu[x]},
#>  {               0}}
```

The mean expectations for the observed variables can be symbolically
derived using the `ramR::g_sym()` function.

``` r
ramR::g_sym(A, u, filter)
#> {{alpha+beta*mu[x]},
#>  {           mu[x]}}
```

## Numerical Example

This is a numerical example for the model

*y* = *α* + *β* ⋅ *x* + *ε*

*y* = 0 + 0.50*x* + *ε*.

``` r
A <- S <- matrix(
  data = 0,
  nrow = 3,
  ncol = 3
)
A[1, ] <- c(0, 1, 1)
diag(S) <- c(0, 0.25, 1)
colnames(A) <- rownames(A) <- c("y", "x", "e")
filter <- diag(2)
filter <- cbind(filter, 0)
colnames(filter) <- c("y", "x", "e")
rownames(filter) <- c("y", "x")
u <- c(0.00, 0.50, 0.00)
```

The covariance expectations can be numerically derived using the
`ramR::C_num()` function.

``` r
ramR::C_num(A, S)
#>      y    x e
#> y 1.25 0.25 1
#> x 0.25 0.25 0
#> e 1.00 0.00 1
```

The covariance expectations for the observed variables can be
numerically derived using the `ramR::M_num()` function.

``` r
ramR::M_num(A, S, filter)
#>      y    x
#> y 1.25 0.25
#> x 0.25 0.25
```

The mean expectations can be numerically derived using the
`ramR::v_num()` function.

``` r
ramR::v_num(A, u)
#>     v
#> y 0.5
#> x 0.5
#> e 0.0
```

The mean expectations for the observed variables can be numerically
derived using the `ramR::v_num()` function.

``` r
ramR::g_num(A, u, filter)
#>     g
#> y 0.5
#> x 0.5
```

## More Information

See [GitHub Pages](https://jeksterslab.github.io/ramR/index.html) for
package documentation. See
[ramR\_notes](https://jeksterslab.github.io/ramR_notes/index.html) for
notes and more examples.
