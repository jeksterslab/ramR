Reticular Action Model (RAM) Notation
================
Ivan Jacob Agaloos Pesigan
2021-01-25

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

## Example

Let *y*, *m*, *x*, *ε*<sub>*y*</sub>, and *ε*<sub>*m*</sub> be random
variables whose associations are given by

*y* = *β*<sub>0</sub> + *β*<sub>1</sub>*x* + *β*<sub>2</sub>*m* + *ε*<sub>*y*</sub>

*m* = *α*<sub>0</sub> + *α*<sub>1</sub>*x* + *ε*<sub>*m*</sub>

where

-   *β*<sub>1</sub> is the path from *x* on *y*
-   *β*<sub>2</sub> is the path from *m* to *y*
-   *α*<sub>1</sub> is the path from *x* to *m*
-   *ε*<sub>*y*</sub> and *ε*<sub>*m*</sub> are uncorrelated error terms
    with means of zero and constant variances of
    *σ*<sub>*ε*<sub>*y*</sub></sub><sup>2</sup> and
    *σ*<sub>*ε*<sub>*m*</sub></sub><sup>2</sup> respectively
-   *β*<sub>0</sub> and *α*<sub>0</sub> are intercepts

### Equations to RAM

The `ramR::eq2ram` function converts structural equations to RAM
notation.

#### Symbolic

The model can be expressed in the following equations

``` r
eq_sym <- "
  # VARIABLE1 OPERATION VARIABLE2 LABEL
  ey          by        y         1
  em          by        m         1
  y           on        x         beta[1]
  y           on        m         beta[2]
  m           on        x         alpha[1]
  ey          with      ey        sigma[varepsilon[y]]^2
  em          with      em        sigma[varepsilon[m]]^2
  x           with      x         sigma[x]^2
  y           on        1         beta[0]
  m           on        1         alpha[0]
  x           on        1         mu[x]
"
```

``` r
ram <- ramR::eq2ram(eq_sym)
ram
#> $eq
#>    var1   op var2                  label
#> 1    ey   by    y                      1
#> 2    em   by    m                      1
#> 3     y   on    x                beta[1]
#> 4     y   on    m                beta[2]
#> 5     m   on    x               alpha[1]
#> 6    ey with   ey sigma[varepsilon[y]]^2
#> 7    em with   em sigma[varepsilon[m]]^2
#> 8     x with    x             sigma[x]^2
#> 9     y   on    1                beta[0]
#> 10    m   on    1               alpha[0]
#> 11    x   on    1                  mu[x]
#> 
#> $variables
#> [1] "y"  "m"  "x"  "ey" "em"
#> 
#> $A
#>    y   m         x          ey  em 
#> y  "0" "beta[2]" "beta[1]"  "1" "0"
#> m  "0" "0"       "alpha[1]" "0" "1"
#> x  "0" "0"       "0"        "0" "0"
#> ey "0" "0"       "0"        "0" "0"
#> em "0" "0"       "0"        "0" "0"
#> 
#> $S
#>    y   m   x            ey                       em                      
#> y  "0" "0" "0"          "0"                      "0"                     
#> m  "0" "0" "0"          "0"                      "0"                     
#> x  "0" "0" "sigma[x]^2" "0"                      "0"                     
#> ey "0" "0" "0"          "sigma[varepsilon[y]]^2" "0"                     
#> em "0" "0" "0"          "0"                      "sigma[varepsilon[m]]^2"
#> 
#> $filter
#>   y m x ey em
#> y 1 0 0  0  0
#> m 0 1 0  0  0
#> x 0 0 1  0  0
#> 
#> $u
#>    u         
#> y  "beta[0]" 
#> m  "alpha[0]"
#> x  "mu[x]"   
#> ey "0"       
#> em "0"
```

#### Numerical

Matrices of the RAM notation can be derived numerically by supplying
values to the parameters.

``` r
eq_num <- "
  # VARIABLE1 OPERATION VARIABLE2 VALUE
  ey          by        y         1
  em          by        m         1
  y           on        x         0.00
  y           on        m         0.50
  m           on        x         0.50
  ey          with      ey        168.75
  em          with      em        168.75
  x           with      x         225
  y           on        1         50
  m           on        1         50
  x           on        1         100
"
```

``` r
ramR::eq2ram(eq_num)
#> $eq
#>    var1   op var2  label
#> 1    ey   by    y   1.00
#> 2    em   by    m   1.00
#> 3     y   on    x   0.00
#> 4     y   on    m   0.50
#> 5     m   on    x   0.50
#> 6    ey with   ey 168.75
#> 7    em with   em 168.75
#> 8     x with    x 225.00
#> 9     y   on    1  50.00
#> 10    m   on    1  50.00
#> 11    x   on    1 100.00
#> 
#> $variables
#> [1] "y"  "m"  "x"  "ey" "em"
#> 
#> $A
#>    y   m   x ey em
#> y  0 0.5 0.0  1  0
#> m  0 0.0 0.5  0  1
#> x  0 0.0 0.0  0  0
#> ey 0 0.0 0.0  0  0
#> em 0 0.0 0.0  0  0
#> 
#> $S
#>    y m   x     ey     em
#> y  0 0   0   0.00   0.00
#> m  0 0   0   0.00   0.00
#> x  0 0 225   0.00   0.00
#> ey 0 0   0 168.75   0.00
#> em 0 0   0   0.00 168.75
#> 
#> $filter
#>   y m x ey em
#> y 1 0 0  0  0
#> m 0 1 0  0  0
#> x 0 0 1  0  0
#> 
#> $u
#>      u
#> y   50
#> m   50
#> x  100
#> ey   0
#> em   0
```

### Equations to Expectations

The `ramR` package has a utility function to convert structural
equations to expectations both symbolically and numerically.

#### Symbolic

``` r
exp <- ramR::eq2exp_sym(eq_sym)
exp
#> $variables
#> [1] "y"  "m"  "x"  "ey" "em"
#> 
#> $A
#> {{       0,  beta[2],  beta[1],        1,        0},
#>  {       0,        0, alpha[1],        0,        1},
#>  {       0,        0,        0,        0,        0},
#>  {       0,        0,        0,        0,        0},
#>  {       0,        0,        0,        0,        0}} 
#> 
#> $S
#> {{                     0,                      0,                      0,                      0,                      0},
#>  {                     0,                      0,                      0,                      0,                      0},
#>  {                     0,                      0,             sigma[x]^2,                      0,                      0},
#>  {                     0,                      0,                      0, sigma[varepsilon[y]]^2,                      0},
#>  {                     0,                      0,                      0,                      0, sigma[varepsilon[m]]^2}} 
#> 
#> $u
#> {{ beta[0]},
#>  {alpha[0]},
#>  {   mu[x]},
#>  {       0},
#>  {       0}} 
#> 
#> $filter
#> {{1, 0, 0, 0, 0},
#>  {0, 1, 0, 0, 0},
#>  {0, 0, 1, 0, 0}} 
#> 
#> $v
#> {{beta[0]+beta[2]*alpha[0]+(beta[1]+beta[2]*alpha[1])*mu[x]},
#>  {                                  alpha[0]+alpha[1]*mu[x]},
#>  {                                                    mu[x]},
#>  {                                                        0},
#>  {                                                        0}} 
#> 
#> $g
#> {{beta[0]+beta[2]*alpha[0]+(beta[1]+beta[2]*alpha[1])*mu[x]},
#>  {                                  alpha[0]+alpha[1]*mu[x]},
#>  {                                                    mu[x]}} 
#> 
#> $C
#> {{sigma[x]^2*(beta[1]+beta[2]*alpha[1])^2+sigma[varepsilon[y]]^2+sigma[varepsilon[m]]^2*beta[2]^2,                   (beta[1]+beta[2]*alpha[1])*sigma[x]^2*alpha[1]+beta[2]*sigma[varepsilon[m]]^2,                                                           (beta[1]+beta[2]*alpha[1])*sigma[x]^2,                                                                          sigma[varepsilon[y]]^2,                                                                  beta[2]*sigma[varepsilon[m]]^2},
#>  {                  alpha[1]*sigma[x]^2*(beta[1]+beta[2]*alpha[1])+sigma[varepsilon[m]]^2*beta[2],                                                    sigma[x]^2*alpha[1]^2+sigma[varepsilon[m]]^2,                                                                             alpha[1]*sigma[x]^2,                                                                                               0,                                                                          sigma[varepsilon[m]]^2},
#>  {                                                          sigma[x]^2*(beta[1]+beta[2]*alpha[1]),                                                                             sigma[x]^2*alpha[1],                                                                                      sigma[x]^2,                                                                                               0,                                                                                               0},
#>  {                                                                         sigma[varepsilon[y]]^2,                                                                                               0,                                                                                               0,                                                                          sigma[varepsilon[y]]^2,                                                                                               0},
#>  {                                                                 sigma[varepsilon[m]]^2*beta[2],                                                                          sigma[varepsilon[m]]^2,                                                                                               0,                                                                                               0,                                                                          sigma[varepsilon[m]]^2}} 
#> 
#> $M
#> {{sigma[x]^2*(beta[1]+beta[2]*alpha[1])^2+sigma[varepsilon[y]]^2+sigma[varepsilon[m]]^2*beta[2]^2,                   (beta[1]+beta[2]*alpha[1])*sigma[x]^2*alpha[1]+beta[2]*sigma[varepsilon[m]]^2,                                                           (beta[1]+beta[2]*alpha[1])*sigma[x]^2},
#>  {                  alpha[1]*sigma[x]^2*(beta[1]+beta[2]*alpha[1])+sigma[varepsilon[m]]^2*beta[2],                                                    sigma[x]^2*alpha[1]^2+sigma[varepsilon[m]]^2,                                                                             alpha[1]*sigma[x]^2},
#>  {                                                          sigma[x]^2*(beta[1]+beta[2]*alpha[1]),                                                                             sigma[x]^2*alpha[1],                                                                                      sigma[x]^2}}
```

#### Numerical

``` r
ramR::eq2exp_num(eq_num)
#> $variables
#> [1] "y"  "m"  "x"  "ey" "em"
#> 
#> $A
#>    y   m   x ey em
#> y  0 0.5 0.0  1  0
#> m  0 0.0 0.5  0  1
#> x  0 0.0 0.0  0  0
#> ey 0 0.0 0.0  0  0
#> em 0 0.0 0.0  0  0
#> 
#> $S
#>    y m   x     ey     em
#> y  0 0   0   0.00   0.00
#> m  0 0   0   0.00   0.00
#> x  0 0 225   0.00   0.00
#> ey 0 0   0 168.75   0.00
#> em 0 0   0   0.00 168.75
#> 
#> $u
#>      u
#> y   50
#> m   50
#> x  100
#> ey   0
#> em   0
#> 
#> $filter
#>   y m x ey em
#> y 1 0 0  0  0
#> m 0 1 0  0  0
#> x 0 0 1  0  0
#> 
#> $v
#>      v
#> y  100
#> m  100
#> x  100
#> ey   0
#> em   0
#> 
#> $g
#>     g
#> y 100
#> m 100
#> x 100
#> 
#> $C
#>          y      m      x     ey      em
#> y  225.000 112.50  56.25 168.75  84.375
#> m  112.500 225.00 112.50   0.00 168.750
#> x   56.250 112.50 225.00   0.00   0.000
#> ey 168.750   0.00   0.00 168.75   0.000
#> em  84.375 168.75   0.00   0.00 168.750
#> 
#> $M
#>        y     m      x
#> y 225.00 112.5  56.25
#> m 112.50 225.0 112.50
#> x  56.25 112.5 225.00
```

## More Information

See [GitHub Pages](https://jeksterslab.github.io/ramR/index.html) for
package documentation. See
[ramR\_notes](https://jeksterslab.github.io/ramR_notes/index.html) for
notes and more examples.
