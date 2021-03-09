#' ---
#' title: "Test: Two-Variable Linear Regression - u"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - u}
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
Numericu <- round(
  ramR::u(A, v),
  digits = 4
)
Symbolicu <- round(
  Ryacas::as_r(
    ramR::u(
      Ryacas::ysym(A),
      v,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicuExpr <- round(
  eval(
    ramR::u(
      Ryacas::ysym(A),
      v,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
#'
#+ str
ramR::u(
  Ryacas::ysym(A),
  v,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  ramR::u(
    Ryacas::ysym(A),
    v,
    R = FALSE,
    format = "tex"
  ),
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
u <- round(u, digits = 4)
#'
#+ testthat
testthat::test_that("u.", {
  for (i in seq_len(nrow(u))) {
    for (j in seq_len(ncol(u))) {
      testthat::expect_equal(
        u[i, j],
        Numericu[i, j],
        Symbolicu[i, j],
        SymbolicuExpr[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ coverage
ramR::u(
  Ryacas::ysym(A),
  v,
  exe = FALSE
)
ramR::u(
  Ryacas::ysym(A),
  Ryacas::ysym(as.vector(v)),
  exe = FALSE
)
ramR::u(
  Ryacas::ysym(A),
  Ryacas::ysym(v),
  exe = FALSE,
  check = FALSE
)
