#' ---
#' title: "Test: Two-Variable Linear Regression - g"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - g}
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
#+
Numericg <- round(
  ramR::g(A, u, Filter),
  digits = 4
)
Numericv <- round(
  ramR::g(A, u, Filter = NULL),
  digits = 4
)
Symbolicg <- round(
  Ryacas::as_r(
    ramR::g(
      Ryacas::ysym(A),
      u,
      Filter,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
Symbolicv <- round(
  Ryacas::as_r(
    ramR::g(
      Ryacas::ysym(A),
      u,
      Filter = NULL,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicgExpr <- round(
  eval(
    ramR::g(
      Ryacas::ysym(A),
      u,
      Filter,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicvExpr <- round(
  eval(
    ramR::g(
      Ryacas::ysym(A),
      u,
      Filter = NULL,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
#'
#+ str
ramR::g(
  Ryacas::ysym(A),
  u,
  Filter,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  ramR::g(
    Ryacas::ysym(A),
    u,
    Filter,
    R = FALSE,
    format = "tex"
  ),
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
g <- round(g, digits = 4)
v <- round(v, digits = 4)
#'
#+ testthat
testthat::test_that("g.", {
  for (i in seq_len(nrow(g))) {
    for (j in seq_len(ncol(g))) {
      testthat::expect_equal(
        g[i, j],
        Numericg[i, j],
        Symbolicg[i, j],
        SymbolicgExpr[i, j],
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
        SymbolicvExpr[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ coverage
ramR::g(
  Ryacas::ysym(A),
  u,
  Ryacas::ysym(Filter),
  exe = FALSE
)
ramR::g(
  Ryacas::ysym(A),
  Ryacas::ysym(as.vector(u)),
  Ryacas::ysym(Filter),
  exe = FALSE
)
ramR::g(
  Ryacas::ysym(A),
  Ryacas::ysym(u),
  Ryacas::ysym(Filter),
  exe = FALSE,
  check = FALSE
)
