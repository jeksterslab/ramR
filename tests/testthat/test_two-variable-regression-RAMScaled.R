#' ---
#' title: "Test: Two-Variable Linear Regression - RAMScaled"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - RAMScaled}
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
NumericRAM <- ramR::RAMScaled(A, S, Filter)
NumericRAM.A.scaled <- NumericRAM$A.scaled
NumericRAM.S.scaled <- NumericRAM$S.scaled
SymbolicRAM <- ramR::RAMScaled(
  Ryacas::ysym(A),
  S,
  Filter,
  R = FALSE,
  simplify = TRUE
)
SymbolicRAM.A.scaled <- Ryacas::as_r(SymbolicRAM$A.scaled)
SymbolicRAM.S.scaled <- Ryacas::as_r(SymbolicRAM$S.scaled)
SymbolicRAMExpr <- ramR::RAMScaled(
  Ryacas::ysym(A),
  S,
  Filter,
  R = TRUE,
  simplify = TRUE
)
SymbolicRAMExpr.A.scaled <- eval(SymbolicRAMExpr$A.scaled)
SymbolicRAMExpr.S.scaled <- eval(SymbolicRAMExpr$S.scaled)
SymbolicRAMExpr2 <- ramR::RAMScaled(
  Ryacas::ysym(A),
  S,
  Filter = NULL,
  R = TRUE,
  simplify = TRUE
)
SymbolicRAMExpr2.A.scaled <- eval(SymbolicRAMExpr2$A.scaled)
SymbolicRAMExpr2.S.scaled <- eval(SymbolicRAMExpr2$S.scaled)
#'
#+ str
ramR::RAMScaled(
  Ryacas::ysym(A),
  S,
  Filter,
  R = FALSE,
  ysym = FALSE,
  simplify = TRUE
)
#'
#+ tex, results = "asis"
latex <- ramR::RAMScaled(
  Ryacas::ysym(A),
  S,
  Filter,
  R = FALSE,
  format = "tex"
)
cat(
  "\\begin{align*}",
  latex$A.scaled,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  latex$S.scaled,
  "\\end{align*}",
  sep = ""
)
#+ round_source
C.scaled <- round(C.scaled, digits = 4)
M.scaled <- round(M.scaled, digits = 4)
#'
#+ testthat
NumericRAM.Expectations <- ramR::Expectations(
  NumericRAM.A.scaled,
  NumericRAM.S.scaled,
  Filter = Filter,
  check = FALSE
)
SymbolicRAM.Expectations <- ramR::Expectations(
  SymbolicRAM.A.scaled,
  SymbolicRAM.S.scaled,
  Filter = Filter,
  check = FALSE
)
SymbolicRAMExpr.Expectations <- ramR::Expectations(
  SymbolicRAMExpr.A.scaled,
  SymbolicRAMExpr.S.scaled,
  Filter = Filter,
  check = FALSE
)
SymbolicRAMExpr2.Expectations <- ramR::Expectations(
  SymbolicRAMExpr2.A.scaled,
  SymbolicRAMExpr2.S.scaled,
  Filter = Filter,
  check = FALSE
)
testthat::test_that("C.scaled", {
  for (i in seq_len(nrow(C.scaled))) {
    for (j in seq_len(ncol(C.scaled))) {
      testthat::expect_equal(
        C.scaled[i, j],
        round(NumericRAM.Expectations$C[i, j], digits = 4),
        round(NumericRAM.Expectations$C.scaled[i, j], digits = 4),
        round(SymbolicRAM.Expectations$C[i, j], digits = 4),
        round(SymbolicRAM.Expectations$C.scaled[i, j], digits = 4),
        round(SymbolicRAMExpr.Expectations$C[i, j], digits = 4),
        round(SymbolicRAMExpr.Expectations$C.scaled[i, j], digits = 4),
        round(SymbolicRAMExpr2.Expectations$C[i, j], digits = 4),
        round(SymbolicRAMExpr2.Expectations$C.scaled[i, j], digits = 4),
        check.attributes = FALSE
      )
    }
  }
})
testthat::test_that("M.scaled", {
  for (i in seq_len(nrow(M.scaled))) {
    for (j in seq_len(ncol(M.scaled))) {
      testthat::expect_equal(
        M.scaled[i, j],
        round(NumericRAM.Expectations$M[i, j], digits = 4),
        round(NumericRAM.Expectations$M.scaled[i, j], digits = 4),
        round(SymbolicRAM.Expectations$M[i, j], digits = 4),
        round(SymbolicRAM.Expectations$M.scaled[i, j], digits = 4),
        round(SymbolicRAMExpr.Expectations$M[i, j], digits = 4),
        round(SymbolicRAMExpr.Expectations$M.scaled[i, j], digits = 4),
        round(SymbolicRAMExpr2.Expectations$M[i, j], digits = 4),
        round(SymbolicRAMExpr2.Expectations$M.scaled[i, j], digits = 4),
        check.attributes = FALSE
      )
    }
  }
})
#'
#+ coverage
ramR::RAMScaled(
  Ryacas::ysym(A),
  Ryacas::ysym(S),
  Ryacas::ysym(Filter),
  exe = FALSE,
  check = FALSE
)
