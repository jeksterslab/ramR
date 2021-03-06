#' ---
#' title: "Test: Two-Variable Linear Regression - S"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - S}
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
NumericS <- round(
  ramR::S(A, C),
  digits = 4
)
SymbolicS <- round(
  Ryacas::as_r(
    ramR::S(
      Ryacas::ysym(A),
      C,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicSExpr <- round(
  eval(
    ramR::S(
      Ryacas::ysym(A),
      C,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
#'
#+ str
ramR::S(
  Ryacas::ysym(A),
  C,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  ramR::S(
    Ryacas::ysym(A),
    C,
    R = FALSE,
    format = "tex"
  ),
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
S <- round(S, digits = 4)
#'
#+ testthat
for (i in seq_len(nrow(S))) {
  for (j in seq_len(ncol(S))) {
    if (abs(NumericS[i, j]) <= 1e-4) {
      NumericS[i, j] <- 0
    }
    if (abs(SymbolicS[i, j]) <= 1e-4) {
      SymbolicS[i, j] <- 0
    }
    if (abs(SymbolicS[i, j]) <= 1e-4) {
      SymbolicS[i, j] <- 0
    }
    if (abs(SymbolicSExpr[i, j]) <= 1e-4) {
      SymbolicSExpr[i, j] <- 0
    }
  }
}
testthat::test_that("S.", {
  for (i in seq_len(nrow(S))) {
    for (j in seq_len(ncol(S))) {
      testthat::expect_equal(
        S[i, j],
        NumericS[i, j],
        SymbolicS[i, j],
        SymbolicSExpr[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ coverage
ramR::S(
  Ryacas::ysym(A),
  Ryacas::ysym(C),
  exe = FALSE,
  check = FALSE
)
