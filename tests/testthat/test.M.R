#' ---
#' title: "Test: M"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: M}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ include = FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#'
#+
library(testthat)
library(ram)
#'
#' ## Parameters
#'
#+ parameters
mux <- mum <- muy <- 100
deltam <- 28.5925808
deltay <- 14.4504452
alpha <- 0.7140742
tauprime <- 0.1414214
beta <- 0.7140742
alphabeta <- alpha * beta
sigma2x <- sigma2m <- sigma2y <- 225
sigma2epsilonm <- 110.2720609
sigma2epsilony <- 73.3220995
parameter_mu <- c(
  x = mux,
  m = mum,
  y = muy
)
parameter_M <- c(
  x = mux,
  m = deltam,
  y = deltay
)
parameter_A <- matrix(
  data = c(
    0,
    alpha,
    tauprime,
    0,
    0,
    beta,
    0,
    0,
    0
  ),
  ncol = 3
)
colnames(parameter_A) <- c("x", "m", "y")
rownames(parameter_A) <- c("x", "m", "y")
parameter_S <- matrix(
  data = c(
    sigma2x,
    0,
    0,
    0,
    sigma2epsilonm,
    0,
    0,
    0,
    sigma2epsilony
  ),
  ncol = 3
)
colnames(parameter_S) <- c("x", "m", "y")
rownames(parameter_S) <- c("x", "m", "y")
parameter_filter <- diag(nrow(parameter_A))
colnames(parameter_filter) <- c("x", "m", "y")
rownames(parameter_filter) <- c("x", "m", "y")
parameter_Sigmatheta <- matrix(
  data = c(
    225,
    160.6667,
    146.5477,
    160.6667,
    225,
    183.3884,
    146.5477,
    183.3884,
    225
  ),
  ncol = 3
)
#'
#' Let $x$, $m$, and $y$ be observed variables
#' whose associations are defined by
#'
#' \begin{equation}
#'   y = `r deltay` + `r tauprime` \cdot x + `r beta` \cdot m + \varepsilon_y
#' \end{equation}
#'
#' and
#'
#' \begin{equation}
#'   m = `r deltam` + `r alpha` \cdot x + \varepsilon_m
#' \end{equation}
#'
#' where
#'
#' \begin{equation}
#'   x \sim \mathcal{N} \left( `r mux`, `r sigma2x` \right) ,
#' \end{equation}
#'
#' \begin{equation}
#'   \varepsilon_y \sim \mathcal{N} \left( 0, \sigma_{\varepsilon_{y}}^{2} = `r sigma2epsilony` \right) ,
#' \end{equation}
#'
#' and
#'
#' \begin{equation}
#'   \varepsilon_m \sim \mathcal{N} \left( 0, \sigma_{\varepsilon_{m}}^{2} = `r sigma2epsilonm` \right) .
#' \end{equation}
#'
#' ### Column Vector of Means and Intercepts
#'
#' \begin{equation}
#'   \mathbf{M}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_M[1]` \\
#'     `r parameter_M[2]` \\
#'     `r parameter_M[3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ### Matrix of Regression Slopes
#'
#' \begin{equation}
#'   \mathbf{A}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_A[1, 1]` & `r parameter_A[1, 2]` & `r parameter_A[1, 3]` \\
#'     `r parameter_A[2, 1]` & `r parameter_A[2, 2]` & `r parameter_A[2, 3]` \\
#'     `r parameter_A[3, 1]` & `r parameter_A[3, 2]` & `r parameter_A[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ### Matrix of Variances
#'
#' \begin{equation}
#'   \mathbf{S}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_S[1, 1]` & `r parameter_S[1, 2]` & `r parameter_S[1, 3]` \\
#'     `r parameter_S[2, 1]` & `r parameter_S[2, 2]` & `r parameter_S[2, 3]` \\
#'     `r parameter_S[3, 1]` & `r parameter_S[3, 2]` & `r parameter_S[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ### Filter Matrix
#'
#' \begin{equation}
#'   \mathbf{F}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_filter[1, 1]` & `r parameter_filter[1, 2]` & `r parameter_filter[1, 3]` \\
#'     `r parameter_filter[2, 1]` & `r parameter_filter[2, 2]` & `r parameter_filter[2, 3]` \\
#'     `r parameter_filter[3, 1]` & `r parameter_filter[3, 2]` & `r parameter_filter[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ### Model-Implied Variance-Covariance Matrix
#'
#' \begin{equation}
#'   \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'   =
#'   \begin{bmatrix}
#'     `r parameter_Sigmatheta[1, 1]` & `r parameter_Sigmatheta[1, 2]` & `r parameter_Sigmatheta[1, 3]` \\
#'     `r parameter_Sigmatheta[2, 1]` & `r parameter_Sigmatheta[2, 2]` & `r parameter_Sigmatheta[2, 3]` \\
#'     `r parameter_Sigmatheta[3, 1]` & `r parameter_Sigmatheta[3, 2]` & `r parameter_Sigmatheta[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ## `mutheta()`
#'
#+
result_M <- M(
  mutheta = parameter_mu,
  A = parameter_A,
  filter = parameter_filter
)
#'
#' ### Column Vector of Means and Intercepts
#'
#' \begin{equation}
#'   \mathbf{M}
#'   =
#'   \begin{bmatrix}
#'     `r result_M[1, 1]` \\
#'     `r result_M[2, 1]` \\
#'     `r result_M[3, 1]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ## testthat
#'
#+
test_that("M.", {
  for (i in seq_along(parameter_M)) {
    expect_equal(
      parameter_M[i],
      result_M[i, 1],
      check.attributes = FALSE,
      tolerance = 0.001
    )
  }
})
