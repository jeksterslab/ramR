#' ---
#' title: "Test: Two-Variable Linear Regression"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression}
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
library(ramR)
library(Ryacas)
#'
#+ parameters
alpha <- runif(n = 1, min = -1, max = 1)
beta <- runif(n = 1, min = -1, max = 1)
sigma2x <- runif(n = 1, min = 0, max = 1)
sigma2epsilon <- runif(n = 1, min = 0, max = 1)
mux <- runif(n = 1, min = -1, max = 1)
muepsilon <- 0
A <- S <- matrix(
  data = 0,
  nrow = 3,
  ncol = 3
)
A[1, 2] <- beta
A[1, 3] <- 1
diag(S) <- c(0, sigma2x, sigma2epsilon)
colnames(A) <- rownames(A) <- c("y", "x", "e")
I <- diag(dim(A)[1])
E <- as_r(solve(ysym(I) - ysym(A)))
filter <- diag(2)
filter <- cbind(filter, 0)
colnames(filter) <- c("y", "x", "e")
rownames(filter) <- c("y", "x")
u <- as.matrix(c(alpha, mux, 0))
v <- as.matrix(c(alpha + (beta * mux), mux, 0.00))
rownames(u) <- rownames(v) <- c("y", "x", "e")
colnames(u) <- "u"
colnames(v) <- "v"
C <- as_r(solve(ysym(I) - ysym(A)) * ysym(S) * t(solve(ysym(I) - ysym(A))))
M <- as_r(ysym(filter) * solve(ysym(I) - ysym(A)) * ysym(S) * t(solve(ysym(I) - ysym(A))) * t(ysym(filter)))
g <- as.matrix(as_r(ysym(filter) * ysym(v)))
#'
#' \begin{equation}
#'   \begin{split}
#'     y &= \alpha + \beta x + \varepsilon \\
#'     y &= `r alpha` + `r beta` x + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+ numeric
results_C_num <- round(C_num(A, S), digits = 4)
results_E_num <- round(E_num(A), digits = 4)
results_M_num <- round(M_num(A, S, filter), digits = 4)
results_S_num <- round(S_num(A, C), digits = 4)
results_v_num <- round(v_num(A, u), digits = 4)
results_u_num <- round(u_num(A, v), digits = 4)
results_g_num <- round(g_num(A, u, filter), digits = 4)
mvn <- mvn(n = 1000, A, S, u, filter, empirical = TRUE)
Sigma.mvn <- round(cov(mvn), digits = 4)
mu.mvn <- as.matrix(round(colMeans(mvn), digits = 4))
model <- "
  # VARIABLE1 OPERATION VARIABLE2 VALUE
  e           by        y         1.00;
  y           on        x         1.00;
  e           with      e         0.25;
  x           with      x         0.25;
  y           on        1         0.00;
  x           on        1         0.50
"
data <- eq2data(model, n = 1000, empirical = TRUE)
Sigma.data <- round(cov(data), digits = 4)
mu.data <- as.matrix(round(colMeans(data), digits = 4))
#'
#+ symbolic
results_C_sym <- round(as_r(C_sym(A, S)), digits = 4)
results_E_sym <- round(as_r(E_sym(A)), digits = 4)
results_M_sym <- round(as_r(M_sym(A, S, filter)), digits = 4)
results_S_sym <- round(as_r(S_sym(A, C)), digits = 4)
results_v_sym <- round(as_r(v_sym(A, u)), digits = 4)
results_u_sym <- round(as_r(u_sym(A, v)), digits = 4)
results_g_sym <- round(as_r(g_sym(A, u, filter)), digits = 4)
#'
#+ round_source
C <- round(C, digits = 4)
E <- round(E, digits = 4)
M <- round(M, digits = 4)
S <- round(S, digits = 4)
v <- round(v, digits = 4)
u <- round(u, digits = 4)
g <- round(g, digits = 4)
#'
#+ testthat
test_that("C.", {
  for (i in seq_len(nrow(C))) {
    for (j in seq_len(ncol(C))) {
      expect_equal(
        C[i, j],
        results_C_num[i, j],
        results_C_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("E.", {
  for (i in seq_len(nrow(E))) {
    for (j in seq_len(ncol(E))) {
      expect_equal(
        E[i, j],
        results_E_num[i, j],
        results_E_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("M.", {
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
      expect_equal(
        M[i, j],
        Sigma.mvn[i, j],
        Sigma.data[i, j],
        results_M_num[i, j],
        results_M_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("S.", {
  for (i in seq_len(nrow(S))) {
    for (j in seq_len(ncol(S))) {
      expect_equal(
        S[i, j],
        results_S_num[i, j],
        results_S_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("v.", {
  for (i in seq_len(nrow(v))) {
    for (j in seq_len(ncol(v))) {
      expect_equal(
        v[i, j],
        results_v_num[i, j],
        results_v_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("u.", {
  for (i in seq_len(nrow(u))) {
    for (j in seq_len(ncol(u))) {
      expect_equal(
        u[i, j],
        results_u_num[i, j],
        results_u_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("g.", {
  for (i in seq_len(nrow(g))) {
    for (j in seq_len(ncol(g))) {
      expect_equal(
        g[i, j],
        mu.mvn[i, j],
        mu.data[i, j],
        results_g_num[i, j],
        results_g_sym[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ eq2ram
A <- S <- matrix(
  data = 0,
  nrow = 3,
  ncol = 3
)
A[1, 2] <- "beta"
A[1, 3] <- 1
diag(S) <- c(0, "sigma[x]^2", "sigma[varepsilon]^2")
filter <- diag(2)
filter <- cbind(filter, 0)
u <- as.matrix(c("alpha", "mu[x]", 0))
model <- "
  # VARIABLE1 OPERATION VARIABLE2 LABEL
  e           by        y         1;
  y           on        x         beta;
  e           with      e         sigma[varepsilon]^2;
  x           with      x         sigma[x]^2;
  y           on        1         alpha;
  x           on        1         mu[x]
"
RAM <- eq2ram(model)
results_A <- RAM$A
results_S <- RAM$S
results_filter <- RAM$filter
results_u <- RAM$u
test_that("A.", {
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(ncol(A))) {
      expect_equal(
        A[i, j],
        results_A[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("S.", {
  for (i in seq_len(nrow(S))) {
    for (j in seq_len(ncol(S))) {
      expect_equal(
        S[i, j],
        results_S[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("filter.", {
  for (i in seq_len(nrow(filter))) {
    for (j in seq_len(ncol(filter))) {
      expect_equal(
        filter[i, j],
        results_filter[i, j],
        check.attributes = FALSE
      )
    }
  }
})
test_that("u.", {
  for (i in seq_len(nrow(u))) {
    for (j in seq_len(ncol(u))) {
      expect_equal(
        u[i, j],
        results_u[i, j],
        check.attributes = FALSE
      )
    }
  }
})
