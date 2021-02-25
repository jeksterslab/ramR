#' ---
#' title: "Test: Two-Variable Linear Regression - EqParse"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - EqParse}
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
#+
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmae2   1.25
    x   with x   sigmax2   0.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
ramR:::EqParse(eq)
#'
#+
eq <- "
  # lhs op   rhs value
    e   by   y   1
    y   on   x   1
    e   with e   1
    x   with x   0.25
    y   on   1   0
    x   on   1   0.50
"
ramR:::EqParse(eq)
#'
#+
# equality constraint
# starting values should be identical
# for parameters constrained to be equal
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
ramR:::EqParse(eq)
#'
#+
# `\n` and `;` as line breaks
eq <- "
    e   by   y   1       NA ; y   on   x   beta    0.00 ;
    e   with e   sigmae2 1.25 ; x   with x   sigmax2 0.25 ;
    y   on   1   alpha   0.00 \n x   on   1   mux     0.50 ;
"
ramR:::EqParse(eq)
#'
#+ malformed_eq
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("incomplete eq - trigger warning in do.call rbind", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
#+ duplicates
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("e by y", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("y on x", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("e with e", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("y on 1", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
#+ arrow_self
eq <- "
  # lhs op   rhs par.label par.start
    e   by   e   1         NA
    e   by   y   1         NA
    y   on   x   beta     0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("e on e", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta1     0.00
    y   on   y   beta2     0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("y on y", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
#+ feedback_loop
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta1     0.00
    x   on   y   beta2     0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
testthat::test_that("y on x and x on y", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
# + duplicates_S
eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta1     0.00
    e   with e   sigmax2   1.25
    x   with x   sigmax2   1.25
    e   with x   sigmaex   1.25
    x   with e   sigmaex   1.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
 "
testthat::test_that("e with x and x with e", {
  testthat::expect_error(
    ramR:::EqParse(eq)
  )
})
