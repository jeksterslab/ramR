#' ---
#' title: "Test: Omega"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Omega}
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
parameter_m <- c(
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
parameter_Omega <- matrix(
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
colnames(parameter_Omega) <- c("x", "m", "y")
rownames(parameter_Omega) <- c("x", "m", "y")
parameter_filter <- diag(nrow(parameter_A))
colnames(parameter_filter) <- c("x", "m", "y")
rownames(parameter_filter) <- c("x", "m", "y")
parameter_Sigma <- matrix(
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
#'   \mathbf{m}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_m[1]` \\
#'     `r parameter_m[2]` \\
#'     `r parameter_m[3]`
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
#'   \boldsymbol{\Omega}
#'   =
#'   \begin{bmatrix}
#'     `r parameter_Omega[1, 1]` & `r parameter_Omega[1, 2]` & `r parameter_Omega[1, 3]` \\
#'     `r parameter_Omega[2, 1]` & `r parameter_Omega[2, 2]` & `r parameter_Omega[2, 3]` \\
#'     `r parameter_Omega[3, 1]` & `r parameter_Omega[3, 2]` & `r parameter_Omega[3, 3]`
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
#'     `r parameter_Sigma[1, 1]` & `r parameter_Sigma[1, 2]` & `r parameter_Sigma[1, 3]` \\
#'     `r parameter_Sigma[2, 1]` & `r parameter_Sigma[2, 2]` & `r parameter_Sigma[2, 3]` \\
#'     `r parameter_Sigma[3, 1]` & `r parameter_Sigma[3, 2]` & `r parameter_Sigma[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ## `mutheta()`
#'
#+
result_Omega <- Omega(
  A = parameter_A,
  Sigmatheta = parameter_Sigma
)
#'
#' ### Matrix of Variances
#'
#' \begin{equation}
#'   \boldsymbol{\Omega}
#'   =
#'   \begin{bmatrix}
#'     `r result_Omega[1, 1]` & `r result_Omega[1, 2]` & `r result_Omega[1, 3]` \\
#'     `r result_Omega[2, 1]` & `r result_Omega[2, 2]` & `r result_Omega[2, 3]` \\
#'     `r result_Omega[3, 1]` & `r result_Omega[3, 2]` & `r result_Omega[3, 3]`
#'   \end{bmatrix}
#' \end{equation}
#'
#' ## testthat
#'
#+
test_that("Omega.", {
  for (i in 1:nrow(parameter_Omega)) {
    for (j in 1:ncol(parameter_Omega)) {
      expect_equal(
        parameter_Omega[i, j],
        result_Omega[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
