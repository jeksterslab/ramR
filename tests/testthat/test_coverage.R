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
#+ E
A <- Omega <- matrix(
  0,
  ncol = 3,
  nrow = 2
)
test_that("E: square.", {
  expect_error(
    E(
      A = A
    )
  )
})
test_that("Sigmatheta: Omega is a square matrix.", {
  expect_error(
    Sigmatheta(
      A = A,
      Omega = Omega
    )
  )
})
A <- matrix(
  0,
  ncol = 3,
  nrow = 3
)
Omega <- matrix(
  0,
  ncol = 2,
  nrow = 2
)
test_that("Sigmatheta: A and Omega have identical dimensions.", {
  expect_error(
    Sigmatheta(
      A = A,
      Omega = Omega
    )
  )
})
A <- matrix(
  data = c(
    0,
    -1,
    -1,
    0
  ),
  nrow = 2,
  ncol = 2
)
test_that("E: inverse.", {
  expect_error(
    E(
      A = A
    )
  )
})
