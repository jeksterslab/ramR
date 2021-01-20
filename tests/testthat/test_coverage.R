#' ---
#' title: "Test: Coverage"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Coverage}
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
#'
#+
library(testthat)
library(ramR)
#'
#+ errors
A <- S <- matrix(
  0,
  ncol = 3,
  nrow = 2
)
test_that("E_num: square.", {
  expect_error(
    E_num(A)
  )
})
test_that("E_sym: square.", {
  expect_error(
    E_sym(A)
  )
})
test_that("C_num: S is a square matrix.", {
  expect_error(
    C_num(A, S)
  )
})
test_that("C_sym: S is a square matrix.", {
  expect_error(
    C_sym(A, S)
  )
})
test_that("S_num: A is a square matrix.", {
  expect_error(
    S_num(A, C)
  )
})
test_that("S_sym: A is a square matrix.", {
  expect_error(
    S_sym(A, C)
  )
})
S <- matrix(
  0,
  ncol = 3,
  nrow = 3
)
test_that("C_num: A and S have the same dimensions.", {
  expect_error(
    C_num(A, S)
  )
})
test_that("C_sym: A and S have the same dimensions.", {
  expect_error(
    C_sym(A, S)
  )
})
A <- matrix(
  0,
  ncol = 3,
  nrow = 3
)
C <- matrix(
  0,
  ncol = 3,
  nrow = 2
)
test_that("S_num: A and C have the same dimensions.", {
  expect_error(
    S_num(A, C)
  )
})
test_that("S_sym: A and C have the same dimensions.", {
  expect_error(
    S_sym(A, C)
  )
})
A <- matrix(
  data = c(0, -1, -1, 0),
  nrow = 2,
  ncol = 2
)
test_that("E_num: singular.", {
  expect_error(
    E_num(A)
  )
})
#'
#+ eq2ram
model <- "
  # VARIABLE1 OPERATION VARIABLE2 LABEL
  e           by        y         1;
  y           on        x         beta;
  e           with      e         sigma[varepsilon]^2;
  x           with      x         sigma[x]^2;
"
eq2ram(model)
