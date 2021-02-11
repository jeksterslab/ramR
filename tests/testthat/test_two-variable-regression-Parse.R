#' ---
#' title: "Test: Two-Variable Linear Regression - Parse"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Two-Variable Linear Regression - Parse}
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
  # lhs op   rhs label   start
    e   by   y   1       NA
    y   on   x   beta    0.00
    e   with e   sigmae2 1.25
    x   with x   sigmax2 0.25
    y   on   1   alpha   0.00
    x   on   1   mux     0.50
"
ramR:::Parse(eq)
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
ramR:::Parse(eq)
#'
#+
# equality constraint
# starting values should be identical
# for parameters constrained to be equal
eq <- "
  # lhs op   rhs label   start
    e   by   y   1       NA
    y   on   x   beta    0.00
    e   with e   sigmax2 1.25
    x   with x   sigmax2 1.25
    y   on   1   alpha   0.00
    x   on   1   mux     0.50
"
ramR:::Parse(eq)
#'
#+
# `\n` and `;` as line breaks
eq <- "
    e   by   y   1       NA ; y   on   x   beta    0.00 ;
    e   with e   sigmae2 1.25 ; x   with x   sigmax2 0.25 ;
    y   on   1   alpha   0.00 \n x   on   1   mux     0.50 ;
"
ramR:::Parse(eq)
