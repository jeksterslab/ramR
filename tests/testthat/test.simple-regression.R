#' ---
#' title: "Test: Simple Regression"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Simple Regression}
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
#' ## Specification 1 - Includes Error Term as a Latent Variable
#'
#+ parameters-01, echo = FALSE
m1 <- -3.951208
m2 <- 13.038328
m3 <- 0
a12 <- 1.269259
a13 <- 1
omega22 <- 7.151261
omega33 <- 47.659854
A <- Omega <- matrix(
  data = 0,
  nrow = 3,
  ncol = 3
)
A[1, 2] <- a12
A[1, 3] <- a13
Omega[2, 2] <- omega22
Omega[3, 3] <- omega33
m <- matrix(
  data = c(
    m1,
    m2,
    m3
  ),
  ncol = 1
)
filter <- diag(3)
muv1 <- m1 + a12 * m2
muv2 <- m2
muv3 <- m3
mu <- matrix(
  data = c(
    muv1,
    muv2,
    muv3
  ),
  ncol = 1
)
cov11 <- a12^2 * omega22 + omega33
cov12 <- a12 * omega22
cov13 <- omega33
cov21 <- cov12
cov22 <- omega22
cov23 <- 0
cov31 <- omega33
cov32 <- cov23
cov33 <- omega33
Cov <- matrix(
  data = c(
    cov11,
    cov21,
    cov31,
    cov12,
    cov22,
    cov32,
    cov13,
    cov23,
    cov33
  ),
  ncol = 3
)
#'
#' Let $v_1$, $v_2$, and $u$ be random variables whose associations are given by the regression equation
#'
#' \begin{equation}
#'   \begin{split}
#'     v_1
#'     &=
#'     m_1 + a_{1, 2} v_2 + u \\
#'     &=
#'     `r m1` + `r a12` \cdot v_2 + u .
#'   \end{split}
#' \end{equation}
#'
#' \noindent $v_1$ and $v_2$ are observed variables and
#' $u$ is a stochastic error term which is normally distributed around zero
#' with constant variance across values of $v_2$
#'
#' \begin{equation}
#'   u
#'   \sim
#'   \mathcal{N} \left( m_3 = 0, \omega_{3, 3} = `r omega33` \right) .
#' \end{equation}
#'
#' \noindent $v_2$ has a mean of $m_2 = `r m2`$ and a variance of $\omega_{2, 2} = `r omega22`$.
#'
#' Below are two ways of specifying this model.
#' The first specification includes the error term $u$ as a latent variable.
#' The second specification only includes the observed variables.
#'
#+ specification-01
m
A
Omega
filter
result_mutheta_01 <- ram::mutheta(
  m,
  A = A,
  filter = filter
)
result_mutheta_01
result_Sigmatheta_01 <- ram::Sigmatheta(
  A = A,
  Omega = Omega,
  filter = filter
)
result_Sigmatheta_01
result_m_01 <- m(
  mutheta = mu,
  A = A,
  filter = filter
)
result_m_01
result_Omega_01 <- Omega(
  A = A,
  Sigmatheta = Cov
)
result_Omega_01
test_that("mutheta-01.", {
  for (i in 1:nrow(mu)) {
    for (j in 1:ncol(mu)) {
      expect_equal(
        mu[i, j],
        result_mutheta_01[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("Sigmatheta-01.", {
  for (i in 1:nrow(Cov)) {
    for (j in 1:ncol(Cov)) {
      expect_equal(
        Cov[i, j],
        result_Sigmatheta_01[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("m-01.", {
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      expect_equal(
        m[i, j],
        result_m_01[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("Omega-01.", {
  for (i in 1:nrow(Omega)) {
    for (j in 1:ncol(Omega)) {
      expect_equal(
        Omega[i, j],
        result_Omega_01[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
#'
#' ## Specification 2 - Observed Variables
#'
#+ parameters-02, echo = FALSE
A <- A[1:2, 1:2]
Omega <- Omega[1:2, 1:2]
Omega[1, 1] <- omega33
omega11 <- Omega[1, 1]
m <- m[1:2, , drop = FALSE]
Cov <- Cov[1:2, 1:2]
mu <- mu[1:2, , drop = FALSE]
filter <- filter[1:2, 1:2]
#'
#+ specification-02
m
A
Omega
result_mutheta_02 <- ram::mutheta(
  m,
  A = A
)
result_mutheta_02
result_Sigmatheta_02 <- ram::Sigmatheta(
  A = A,
  Omega = Omega
)
result_Sigmatheta_02
result_m_02 <- m(
  mutheta = mu,
  A = A
)
result_m_02
result_Omega_02 <- Omega(
  A = A,
  Sigmatheta = Cov
)
result_Omega_02
test_that("mutheta-02.", {
  for (i in 1:nrow(mu)) {
    for (j in 1:ncol(mu)) {
      expect_equal(
        mu[i, j],
        result_mutheta_02[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("Sigmatheta-02.", {
  for (i in 1:nrow(Cov)) {
    for (j in 1:ncol(Cov)) {
      expect_equal(
        Cov[i, j],
        result_Sigmatheta_02[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("m-02.", {
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      expect_equal(
        m[i, j],
        result_m_02[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
test_that("Omega-02.", {
  for (i in 1:nrow(Omega)) {
    for (j in 1:ncol(Omega)) {
      expect_equal(
        Omega[i, j],
        result_Omega_02[i, j],
        check.attributes = FALSE,
        tolerance = 0.001
      )
    }
  }
})
#'
#+ coverage
m <- as.vector(m)
result_mutheta_03 <- as.vector(
  ram::mutheta(
    m = m,
    A = A
  )
)
mu <- as.vector(mu)
result_m_03 <- as.vector(
  ram::m(
    mutheta = mu,
    A = A
  )
)
test_that("mutheta-03.", {
  for (i in seq_along(mu)) {
    expect_equal(
      mu[i],
      result_mutheta_03[i],
      check.attributes = FALSE,
      tolerance = 0.001
    )
  }
})
test_that("m-03.", {
  for (i in seq_along(m)) {
    expect_equal(
      m[i],
      result_m_03[i],
      check.attributes = FALSE,
      tolerance = 0.001
    )
  }
})
