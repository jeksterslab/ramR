#' ---
#' title: "Test: Two-Variable Linear Regression - Eq2Expectations"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - Eq2Expectations}
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
#' ## ASCII Characters
#'
#+ ascii
eq <- "
  # lhs op   rhs par.label
    e   by   y   1
    y   on   x   beta
    e   with e   sigmae2
    x   with x   sigmax2
    y   on   1   alpha
    x   on   1   mux
"
SymbolicLabels <- ramR::Eq2Expectations(eq, par = FALSE, R = TRUE)
SymbolicPars <- ramR::Eq2Expectations(eq, par = TRUE, R = TRUE)
#'
#' ## LaTeX
#'
#+ latex
eq <- "
  # lhs op   rhs par.label
    e   by   y   1
    y   on   x   beta
    e   with e   sigma[epsilon]^2
    x   with x   sigma[x]^2
    y   on   1   alpha
    x   on   1   mu[x]
"
SymbolicLabelsTex <- ramR::Eq2Expectations(
  eq,
  par = FALSE,
  R = FALSE,
  format = "tex"
)
SymbolicParsTex <- ramR::Eq2Expectations(
  eq,
  par = TRUE,
  R = FALSE,
  format = "tex"
)
#'
#+ SymbolicLabelsTex, results = "asis"
cat(
  "\\begin{align*}",
  "\\mathbf{v} = ",
  SymbolicLabelsTex$v,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{g} = ",
  SymbolicLabelsTex$g,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{C} = ",
  SymbolicLabelsTex$C,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{M} = ",
  SymbolicLabelsTex$M,
  "\\end{align*}",
  sep = ""
)
#'
#+ SymbolicParsTex, results = "asis"
cat(
  "\\begin{align*}",
  "\\mathbf{v} = ",
  SymbolicParsTex$v,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{g} = ",
  SymbolicParsTex$g,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{C} = ",
  SymbolicParsTex$C,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{M} = ",
  SymbolicParsTex$M,
  "\\end{align*}",
  sep = ""
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
#'     y &= `r alpha` + \left( `r beta` \right) x + \varepsilon
#'   \end{split}
#' \end{equation}
#'
#+ SymbolicLabels
SymbolicLabelsv <- round(eval(SymbolicLabels$v), digits = 4)
SymbolicLabelsg <- round(eval(SymbolicLabels$g), digits = 4)
SymbolicLabelsC <- round(eval(SymbolicLabels$C), digits = 4)
SymbolicLabelsM <- round(eval(SymbolicLabels$M), digits = 4)
#'
#+ SymbolicPars
p <- c(beta, sigmae2, sigmax2, alpha, mux)
p1 <- p[1]
p2 <- p[2]
p3 <- p[3]
p4 <- p[4]
p5 <- p[5]
SymbolicParsv <- round(eval(SymbolicPars$v), digits = 4)
SymbolicParsg <- round(eval(SymbolicPars$g), digits = 4)
SymbolicParsC <- round(eval(SymbolicPars$C), digits = 4)
SymbolicParsM <- round(eval(SymbolicPars$M), digits = 4)
#'
#+ numeric
eq <- paste(
  "e by y 1", "\n",
  "y on x", beta, "\n",
  "e with e", sigmae2, "\n",
  "x with x", sigmax2, "\n",
  "y on 1", alpha, "\n",
  "x on 1", mux, "\n"
)
NumericLabels <- ramR::Eq2Expectations(eq, par = FALSE)
NumericPars <- ramR::Eq2Expectations(eq, par = TRUE)
#'
#+ NumericLabels
NumericLabelsv <- round(NumericLabels$v, digits = 4)
NumericLabelsg <- round(NumericLabels$g, digits = 4)
NumericLabelsC <- round(NumericLabels$C, digits = 4)
NumericLabelsM <- round(NumericLabels$M, digits = 4)
#'
#+ NumericPars
NumericParsv <- round(NumericPars$v, digits = 4)
NumericParsg <- round(NumericPars$g, digits = 4)
NumericParsC <- round(NumericPars$C, digits = 4)
NumericParsM <- round(NumericPars$M, digits = 4)
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
        NumericLabelsC[i, j],
        NumericParsC[i, j],
        SymbolicLabelsC[i, j],
        SymbolicParsC[i, j],
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
        NumericLabelsM[i, j],
        NumericParsM[i, j],
        SymbolicLabelsM[i, j],
        SymbolicParsM[i, j],
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
        NumericLabelsv[i, j],
        NumericParsv[i, j],
        SymbolicLabelsv[i, j],
        SymbolicParsv[i, j],
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
        NumericLabelsg[i, j],
        NumericParsg[i, j],
        SymbolicLabelsg[i, j],
        SymbolicParsg[i, j],
        check.attributes = FALSE
      )
    }
  }
})
