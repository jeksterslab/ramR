#' ---
#' title: "Test: Two-Variable Linear Regression - M"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - M}
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
#'     y &= `r alpha` + \left( `r beta` x \right) + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+
NumericM <- round(
  ramR::M(A, S, Filter),
  digits = 4
)
NumericC <- round(
  ramR::M(A, S, Filter = NULL),
  digits = 4
)
SymbolicM <- round(
  Ryacas::as_r(
    ramR::M(
      Ryacas::ysym(A),
      S,
      Filter,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicC <- round(
  Ryacas::as_r(
    ramR::M(
      Ryacas::ysym(A),
      S,
      Filter = NULL,
      R = FALSE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicMExpr <- round(
  eval(
    ramR::M(
      Ryacas::ysym(A),
      S,
      Filter,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
SymbolicCExpr <- round(
  eval(
    ramR::M(
      Ryacas::ysym(A),
      S,
      Filter = NULL,
      R = TRUE,
      simplify = TRUE
    )
  ),
  digits = 4
)
#'
#+ str
ramR::M(
  Ryacas::ysym(A),
  S,
  Filter,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
cat(
  "\\begin{align*}",
  ramR::M(
    Ryacas::ysym(A),
    S,
    Filter,
    R = FALSE,
    format = "tex"
  ),
  "\\end{align*}",
  sep = ""
)
#'
#+ round_source
M <- round(M, digits = 4)
C <- round(C, digits = 4)
#'
#+ testthat
testthat::test_that("M.", {
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
      testthat::expect_equal(
        M[i, j],
        NumericM[i, j],
        SymbolicM[i, j],
        SymbolicMExpr[i, j],
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
        SymbolicCExpr[i, j],
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ coverage
ramR::M(
  Ryacas::ysym(A),
  Ryacas::ysym(S),
  Ryacas::ysym(Filter),
  exe = FALSE
)
