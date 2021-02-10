#' ---
#' title: "Test: Two-Variable Linear Regression - Expectations"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - Expectations}
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
#+ parameters
alpha <- runif(n = 1, min = -1, max = 1)
beta <- runif(n = 1, min = -1, max = 1)
sigmax2 <- runif(n = 1, min = 0, max = 1)
sigmae2 <- runif(n = 1, min = 0, max = 1)
mux <- runif(n = 1, min = -1, max = 1)
muepsilon <- 0
A <- S <- matrixR::ZeroMatrix(3)
A[1, ] <- c(0, beta, 1)
diag(S) <- c(0, sigmax2, sigmae2)
colnames(A) <- rownames(A) <- c("y", "x", "e")
I <- matrixR::IdentityFrom(A)
IminusA <- I - A
E <- Ryacas::as_r(solve(Ryacas::ysym(I) - Ryacas::ysym(A)))
Filter <- diag(2)
Filter <- cbind(Filter, 0)
colnames(Filter) <- c("y", "x", "e")
rownames(Filter) <- c("y", "x")
u <- as.matrix(c(alpha, mux, 0))
v <- as.matrix(c(alpha + (beta * mux), mux, 0.00))
rownames(u) <- rownames(v) <- c("y", "x", "e")
colnames(u) <- "u"
colnames(v) <- "v"
C <- Ryacas::as_r(solve(Ryacas::ysym(I) - Ryacas::ysym(A)) * Ryacas::ysym(S) * t(solve(Ryacas::ysym(I) - Ryacas::ysym(A))))
M <- Ryacas::as_r(Ryacas::ysym(Filter) * solve(Ryacas::ysym(I) - Ryacas::ysym(A)) * Ryacas::ysym(S) * t(solve(Ryacas::ysym(I) - Ryacas::ysym(A))) * t(Ryacas::ysym(Filter)))
g <- as.matrix(Ryacas::as_r(Ryacas::ysym(Filter) * Ryacas::ysym(v)))
#'
#' \begin{equation}
#'   \begin{split}
#'     y &= \alpha + \beta x + \varepsilon \\
#'     y &= `r alpha` + `r beta` x + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+ numeric
Numeric <- ramR::Expectations(A, S, u, Filter)
NumericA <- round(Numeric$A, digits = 4)
NumericS <- round(Numeric$S, digits = 4)
Numericu <- round(Numeric$u, digits = 4)
NumericFilter <- round(Numeric$Filter, digits = 4)
NumericC <- round(Numeric$C, digits = 4)
NumericM <- round(Numeric$M, digits = 4)
Numericv <- round(Numeric$v, digits = 4)
Numericg <- round(Numeric$g, digits = 4)
#'
#+ symbolic
Symbolic <- ramR::Expectations(Ryacas::ysym(A), S, u, Filter)
SymbolicA <- round(Ryacas::as_r(Symbolic$A), digits = 4)
SymbolicS <- round(Ryacas::as_r(Symbolic$S), digits = 4)
Symbolicu <- round(Ryacas::as_r(Symbolic$u), digits = 4)
SymbolicFilter <- round(Ryacas::as_r(Symbolic$Filter), digits = 4)
SymbolicC <- round(Ryacas::as_r(Symbolic$C), digits = 4)
SymbolicM <- round(Ryacas::as_r(Symbolic$M), digits = 4)
Symbolicv <- round(Ryacas::as_r(Symbolic$v), digits = 4)
Symbolicg <- round(Ryacas::as_r(Symbolic$g), digits = 4)
SymbolicExpr <- ramR::Expectations(Ryacas::ysym(A), S, u, Filter, str = FALSE, simplify = TRUE)
SymbolicExprC <- round(Ryacas::as_r(SymbolicExpr$C), digits = 4)
SymbolicExprM <- round(Ryacas::as_r(SymbolicExpr$M), digits = 4)
SymbolicExprv <- round(Ryacas::as_r(SymbolicExpr$v), digits = 4)
SymbolicExprg <- round(Ryacas::as_r(SymbolicExpr$g), digits = 4)
ramR::Expectations(Ryacas::ysym(A), S, u, Filter, str = TRUE, ysym = FALSE, simplify = TRUE)
ExpectationsTex <- ramR::Expectations(Ryacas::ysym(A), S, u, Filter, str = TRUE, tex = TRUE, simplify = TRUE)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  "\\mathbf{A} = ",
  ExpectationsTex$A,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{S} = ",
  ExpectationsTex$S,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{u} = ",
  ExpectationsTex$u,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{F} = ",
  ExpectationsTex$Filter,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{v} = ",
  ExpectationsTex$v,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{g} = ",
  ExpectationsTex$g,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{C} = ",
  ExpectationsTex$C,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{M} = ",
  ExpectationsTex$M,
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
A <- round(A, digits = 4)
S <- round(S, digits = 4)
u <- round(u, digits = 4)
Filter <- round(Filter, digits = 4)
C <- round(C, digits = 4)
M <- round(M, digits = 4)
v <- round(v, digits = 4)
g <- round(g, digits = 4)
#'
#+ testthat
testthat::test_that("A.", {
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(ncol(A))) {
      testthat::expect_equal(
        A[i, j],
        NumericA[i, j],
        SymbolicA[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("S.", {
  for (i in seq_len(nrow(S))) {
    for (j in seq_len(ncol(S))) {
      testthat::expect_equal(
        S[i, j],
        NumericS[i, j],
        SymbolicS[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("u.", {
  for (i in seq_len(nrow(u))) {
    for (j in seq_len(ncol(u))) {
      testthat::expect_equal(
        u[i, j],
        Numericu[i, j],
        Symbolicu[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("Filter.", {
  for (i in seq_len(nrow(Filter))) {
    for (j in seq_len(ncol(Filter))) {
      testthat::expect_equal(
        Filter[i, j],
        NumericFilter[i, j],
        SymbolicFilter[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("C.", {
  for (i in seq_len(nrow(C))) {
    for (j in seq_len(ncol(C))) {
      testthat::expect_equal(
        C[i, j],
        NumericC[i, j],
        SymbolicC[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("M.", {
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
      testthat::expect_equal(
        M[i, j],
        NumericM[i, j],
        SymbolicM[i, j],
        SymbolicExprM[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("v.", {
  for (i in seq_len(nrow(v))) {
    for (j in seq_len(ncol(v))) {
      testthat::expect_equal(
        v[i, j],
        Numericv[i, j],
        Symbolicv[i, j],
        SymbolicExprv[i, j],
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("g.", {
  for (i in seq_len(nrow(g))) {
    for (j in seq_len(ncol(g))) {
      testthat::expect_equal(
        g[i, j],
        Numericg[i, j],
        Symbolicg[i, j],
        SymbolicExprg[i, j],
        check.attributes = FALSE
      )
    }
  }
})
