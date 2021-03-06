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
C <- Ryacas::as_r(
  solve(
    Ryacas::ysym(I) - Ryacas::ysym(A)
  ) * Ryacas::ysym(S) * t(
    solve(Ryacas::ysym(I) - Ryacas::ysym(A))
  )
)
C.scaled <- stats::cov2cor(C)
M <- Ryacas::as_r(
  Ryacas::ysym(Filter) * solve(
    Ryacas::ysym(I) - Ryacas::ysym(A)
  ) * Ryacas::ysym(S) * t(
    solve(Ryacas::ysym(I) - Ryacas::ysym(A))
  ) * t(Ryacas::ysym(Filter))
)
M.scaled <- Filter %*% C.scaled %*% t(Filter)
g <- as.matrix(
  Ryacas::as_r(
    Ryacas::ysym(Filter) * Ryacas::ysym(v)
  )
)
#'
#' \begin{equation}
#'   \begin{split}
#'     y &= \alpha + \beta x + \varepsilon \\
#'     y &= `r alpha` + \left( `r beta` x \right) + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+ numeric
Numeric <- ramR::Expectations(A, S, u, Filter)
NumericC <- round(Numeric$C, digits = 4)
NumericM <- round(Numeric$M, digits = 4)
Numericv <- round(Numeric$v, digits = 4)
Numericg <- round(Numeric$g, digits = 4)
#'
#+ numeric_u_NULL
NumericuNull <- ramR::Expectations(A, S, u = NULL, Filter)
NumericuNullC <- round(NumericuNull$C, digits = 4)
NumericuNullM <- round(NumericuNull$M, digits = 4)
NumericuNullv <- NumericuNull$v
NumericuNullg <- NumericuNull$g
#'
#+ numeric_Filter_NULL
NumericFilterNull <- ramR::Expectations(A, S, u, Filter = NULL)
NumericFilterNullC <- round(NumericFilterNull$C, digits = 4)
NumericFilterNullM <- round(NumericFilterNull$M, digits = 4)
NumericFilterNullv <- round(NumericFilterNull$v, digits = 4)
NumericFilterNullg <- round(NumericFilterNull$g, digits = 4)
#'
#+ symbolic
Symbolic <- ramR::Expectations(Ryacas::ysym(A), S, u, Filter)
SymbolicC <- round(Ryacas::as_r(Symbolic$C), digits = 4)
SymbolicM <- round(Ryacas::as_r(Symbolic$M), digits = 4)
Symbolicv <- round(Ryacas::as_r(Symbolic$v), digits = 4)
Symbolicg <- round(Ryacas::as_r(Symbolic$g), digits = 4)
SymbolicExpr <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u,
  Filter,
  R = TRUE,
  simplify = TRUE
)
SymbolicExprC <- round(Ryacas::as_r(SymbolicExpr$C), digits = 4)
SymbolicExprM <- round(Ryacas::as_r(SymbolicExpr$M), digits = 4)
SymbolicExprv <- round(Ryacas::as_r(SymbolicExpr$v), digits = 4)
SymbolicExprg <- round(Ryacas::as_r(SymbolicExpr$g), digits = 4)
ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u,
  Filter,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u = NULL,
  Filter = NULL,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u = NULL,
  Filter = NULL,
  R = FALSE,
  format = "tex",
  simplify = TRUE
)
ExpectationsTex <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u,
  Filter,
  R = FALSE,
  format = "tex",
  simplify = TRUE
)
#'
#+ symbolic_u_null
SymbolicuNull <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u = NULL,
  Filter
)
SymbolicuNullC <- round(Ryacas::as_r(SymbolicuNull$C), digits = 4)
SymbolicuNullM <- round(Ryacas::as_r(SymbolicuNull$M), digits = 4)
SymbolicuNullv <- SymbolicuNull$v
SymbolicuNullg <- SymbolicuNull$g
SymbolicuNullExpr <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u = NULL,
  Filter,
  R = TRUE,
  simplify = TRUE
)
SymbolicuNullExprC <- round(
  Ryacas::as_r(
    SymbolicuNullExpr$C
  ),
  digits = 4
)
SymbolicuNullExprM <- round(Ryacas::as_r(SymbolicExpr$M), digits = 4)
SymbolicuNullExprv <- SymbolicuNullExpr$v
SymbolicuNullExprg <- SymbolicuNullExpr$g
#'
#+ symbolic_Filter_NULL
SymbolicFilterNull <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u,
  Filter = NULL
)
SymbolicFilterNullC <- round(
  Ryacas::as_r(
    SymbolicFilterNull$C
  ),
  digits = 4
)
SymbolicFilterNullM <- round(
  Ryacas::as_r(
    SymbolicFilterNull$M
  ),
  digits = 4
)
SymbolicFilterNullv <- round(
  Ryacas::as_r(
    SymbolicFilterNull$v
  ),
  digits = 4
)
SymbolicFilterNullg <- round(
  Ryacas::as_r(
    SymbolicFilterNull$g
  ),
  digits = 4
)
SymbolicFilterNullExpr <- ramR::Expectations(
  Ryacas::ysym(A),
  S,
  u,
  Filter = NULL,
  R = TRUE,
  simplify = TRUE
)
SymbolicFilterNullExprC <- round(
  Ryacas::as_r(
    SymbolicFilterNullExpr$C
  ),
  digits = 4
)
SymbolicFilterNullExprM <- round(
  Ryacas::as_r(
    SymbolicFilterNullExpr$M
  ),
  digits = 4
)
SymbolicFilterNullExprv <- round(
  Ryacas::as_r(
    SymbolicFilterNullExpr$v
  ),
  digits = 4
)
SymbolicFilterNullExprg <- round(
  Ryacas::as_r(
    SymbolicFilterNullExpr$g
  ),
  digits = 4
)
#'
#+ tex, results = "asis"
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
C <- round(C, digits = 4)
M <- round(M, digits = 4)
v <- round(v, digits = 4)
g <- round(g, digits = 4)
#'
#+ testthat
testthat::test_that("C.", {
  for (i in seq_len(nrow(C))) {
    for (j in seq_len(ncol(C))) {
      testthat::expect_equal(
        C[i, j],
        NumericC[i, j],
        NumericuNullC[i, j],
        NumericFilterNullC[i, j],
        NumericFilterNullM[i, j],
        SymbolicC[i, j],
        SymbolicuNullC[i, j],
        SymbolicExprC[i, j],
        SymbolicuNullExprC[i, j],
        SymbolicFilterNullC[i, j],
        SymbolicFilterNullM[i, j],
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
        NumericuNullM[i, j],
        SymbolicM[i, j],
        SymbolicuNullM[i, j],
        SymbolicExprM[i, j],
        SymbolicuNullExprM[i, j],
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
        NumericFilterNullv[i, j],
        Symbolicv[i, j],
        SymbolicExprv[i, j],
        SymbolicFilterNullv[i, j],
        SymbolicFilterNullExprv[i, j],
        check.attributes = FALSE
      )
    }
  }
  testthat::expect_true(
    is.null(NumericuNullv)
  )
  testthat::expect_true(
    is.null(SymbolicuNullv)
  )
  testthat::expect_true(
    is.null(SymbolicuNullExprv)
  )
})
testthat::test_that("g.", {
  for (i in seq_len(nrow(g))) {
    for (j in seq_len(ncol(g))) {
      testthat::expect_equal(
        g[i, j],
        Numericg[i, j],
        NumericFilterNullg[i, j],
        Symbolicg[i, j],
        SymbolicExprg[i, j],
        SymbolicFilterNullg[i, j],
        SymbolicFilterNullExprg[i, j],
        check.attributes = FALSE
      )
    }
  }
  testthat::expect_true(
    is.null(NumericuNullg)
  )
  testthat::expect_true(
    is.null(SymbolicuNullg)
  )
  testthat::expect_true(
    is.null(SymbolicuNullExprg)
  )
})
#'
#+ coverage
Symbolic <- ramR::Expectations(
  Ryacas::ysym(A),
  Ryacas::ysym(S),
  Ryacas::ysym(as.vector(u)),
  Ryacas::ysym(Filter)
)
Symbolic <- ramR::Expectations(
  Ryacas::ysym(A),
  Ryacas::ysym(S),
  Ryacas::ysym(u),
  Ryacas::ysym(Filter),
  check = FALSE
)
