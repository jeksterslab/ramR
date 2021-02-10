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
  # lhs op   rhs label
    e   by   y   1
    y   on   x   beta
    e   with e   sigmae2
    x   with x   sigmax2
    y   on   1   alpha
    x   on   1   mux
"
SymbolicLabels <- ramR::Eq2Expectations(eq, par = FALSE, str = FALSE)
SymbolicPars <- ramR::Eq2Expectations(eq, par = TRUE, str = FALSE)
#'
#' ## LaTeX
#'
#+ latex
eq <- "
  # lhs op   rhs label
    e   by   y   1
    y   on   x   beta
    e   with e   sigma[epsilon]^2
    x   with x   sigma[x]^2
    y   on   1   alpha
    x   on   1   mu[x]
"
SymbolicLabelsTex <- ramR::Eq2Expectations(eq, par = FALSE, str = TRUE, tex = TRUE)
SymbolicParsTex <- ramR::Eq2Expectations(eq, par = TRUE, str = TRUE, tex = TRUE)
#'
#+ SymbolicLabelsTex, results = "asis"
cat(
  "\\begin{align*}",
  "\\mathbf{A} = ",
  SymbolicLabelsTex$A,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{S} = ",
  SymbolicLabelsTex$S,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{u} = ",
  SymbolicLabelsTex$u,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{F} = ",
  SymbolicLabelsTex$Filter,
  "\\end{align*}",
  sep = ""
)
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
  "\\mathbf{A} = ",
  SymbolicParsTex$A,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{S} = ",
  SymbolicParsTex$S,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{u} = ",
  SymbolicParsTex$u,
  "\\end{align*}",
  sep = ""
)
cat(
  "\\begin{align*}",
  "\\mathbf{F} = ",
  SymbolicParsTex$Filter,
  "\\end{align*}",
  sep = ""
)
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
#+ SymbolicLabels
SymbolicLabelsA <- round(eval(SymbolicLabels$A), digits = 4)
SymbolicLabelsS <- round(eval(SymbolicLabels$S), digits = 4)
SymbolicLabelsu <- round(eval(SymbolicLabels$u), digits = 4)
SymbolicLabelsFilter <- round(eval(SymbolicLabels$Filter), digits = 4)
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
SymbolicParsA <- round(eval(SymbolicPars$A), digits = 4)
SymbolicParsS <- round(eval(SymbolicPars$S), digits = 4)
SymbolicParsu <- round(eval(SymbolicPars$u), digits = 4)
SymbolicParsFilter <- round(eval(SymbolicPars$Filter), digits = 4)
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
NumericLabelsA <- round(NumericLabels$A, digits = 4)
NumericLabelsS <- round(NumericLabels$S, digits = 4)
NumericLabelsu <- round(NumericLabels$u, digits = 4)
NumericLabelsFilter <- round(NumericLabels$Filter, digits = 4)
NumericLabelsv <- round(NumericLabels$v, digits = 4)
NumericLabelsg <- round(NumericLabels$g, digits = 4)
NumericLabelsC <- round(NumericLabels$C, digits = 4)
NumericLabelsM <- round(NumericLabels$M, digits = 4)
#'
#+ NumericPars
NumericParsA <- round(NumericPars$A, digits = 4)
NumericParsS <- round(NumericPars$S, digits = 4)
NumericParsu <- round(NumericPars$u, digits = 4)
NumericParsFilter <- round(NumericPars$Filter, digits = 4)
NumericParsv <- round(NumericPars$v, digits = 4)
NumericParsg <- round(NumericPars$g, digits = 4)
NumericParsC <- round(NumericPars$C, digits = 4)
NumericParsM <- round(NumericPars$M, digits = 4)
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
        NumericLabelsA[i, j],
        NumericParsA[i, j],
        SymbolicLabelsA[i, j],
        SymbolicParsA[i, j],
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
        NumericLabelsS[i, j],
        NumericParsS[i, j],
        SymbolicLabelsS[i, j],
        SymbolicParsS[i, j],
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
        NumericLabelsu[i, j],
        NumericParsu[i, j],
        SymbolicLabelsu[i, j],
        SymbolicParsu[i, j],
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
        NumericLabelsFilter[i, j],
        NumericParsFilter[i, j],
        SymbolicLabelsFilter[i, j],
        SymbolicParsFilter[i, j],
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
