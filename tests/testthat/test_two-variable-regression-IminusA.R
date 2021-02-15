#' ---
#' title: "Test: Two-Variable Linear Regression - IminusA"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - IminusA}
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
M <- Ryacas::as_r(
  Ryacas::ysym(Filter) * solve(
    Ryacas::ysym(I) - Ryacas::ysym(A)
  ) * Ryacas::ysym(S) * t(
    solve(Ryacas::ysym(I) - Ryacas::ysym(A))
  ) * t(Ryacas::ysym(Filter))
)
g <- as.matrix(
  Ryacas::as_r(
    Ryacas::ysym(Filter) * Ryacas::ysym(v)
  )
)
#'
#' \begin{equation}
#'   \begin{split}
#'     y &= \alpha + \beta x + \varepsilon \\
#'     y &= `r alpha` + `r beta` x + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+
NumericIminusA <- round(
  ramR::IminusA(A),
  digits = 4
)
SymbolicIminusA <- round(
  Ryacas::as_r(
    ramR::IminusA(
      Ryacas::ysym(A),
      str = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicIminusAExpr <- round(
  eval(
    ramR::IminusA(
      Ryacas::ysym(A),
      str = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
#'
#+ str
ramR::IminusA(
  Ryacas::ysym(A),
  str = TRUE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  ramR::IminusA(
    Ryacas::ysym(A),
    str = TRUE,
    tex = TRUE
  ),
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
IminusA <- round(IminusA, digits = 4)
#'
#+ testthat
testthat::test_that("IminusA.", {
  for (i in seq_len(nrow(IminusA))) {
    for (j in seq_len(ncol(IminusA))) {
      testthat::expect_equal(
        IminusA[i, j],
        NumericIminusA[i, j],
        SymbolicIminusA[i, j],
        SymbolicIminusAExpr[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+
ramR::IminusA(
  Ryacas::ysym(A),
  exe = FALSE
)
