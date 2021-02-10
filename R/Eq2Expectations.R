#' Equations to Expectations
#'
#' Converts equations to expectations.
#'
#' The input is a character string
#'   that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs operation rhs label`
#'
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{`left-hand side` measured **by** `right-hand side`}
#'     \item{on}{`left-hand side` regressed **on** `right-hand side`}
#'     \item{with}{`left-hand side` covarying **with** `right-hand side`}
#'     \item{on 1}{`left-hand side` regressed **on 1** for mean structure}
#'   }
#'
#' @section label:
#'   Each parameter should be labeled.
#'   The `label` should be a number for fixed parameters
#'   and a character string for free parameters.
#'   Equality contraints can be imposed by using the same label.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inherit ramR references
#' @inheritParams Eq2RAM
#' @inheritParams Expectations.yac_symbol
#' @examples
#' # Numerical -------------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs label
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' 
#' # Symbolic -------------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs label
#'     e   by   y   1
#'     y   on   x   beta
#'     e   with e   sigmae2
#'     x   with x   sigmax2
#'     y   on   1   alpha
#'     x   on   1   mux
#' "
#' Eq2Expectations(eq, par = FALSE)
#'
#' # Expressions using `label`
#'
#' beta <- 1
#' sigmae2 <- 1
#' sigmax2 <- 0.25
#' alpha <- 0
#' mux <- 0.50
#' Exp <- Eq2Expectations(eq, par = FALSE, str = FALSE)
#' eval(Exp$M)
#' eval(Exp$g)
#'
#' # Expressions using `par.index`
#'
#' p <- c(beta, sigmae2, sigmax2, alpha, mux)
#' p1 <- p[1]
#' p2 <- p[2]
#' p3 <- p[3]
#' p4 <- p[4]
#' p5 <- p[5]
#' Exp <- Eq2Expectations(eq, par = TRUE, str = FALSE)
#' eval(Exp$M)
#' eval(Exp$g)
#' @export
Eq2Expectations <- function(eq,
                            par = FALSE,
                            str = TRUE,
                            ysym = TRUE,
                            simplify = FALSE,
                            tex = FALSE) {
  Expectations <- Eq2RAM(eq, par = par)
  variables <- list(
    variables = Expectations[["variables"]]
  )
  if (is.numeric(Expectations$eq$label)) {
    out <- Expectations.default(
      A = Expectations$A,
      S = Expectations$S,
      u = Expectations$u,
      Filter = Expectations$Filter
    )
  } else {
    out <- Expectations.yac_symbol(
      A = Ryacas::ysym(Expectations$A),
      S = Expectations$S,
      u = Expectations$u,
      Filter = Expectations$Filter,
      str = str,
      ysym = ysym,
      simplify = simplify,
      tex = tex
    )
  }
  return(
    c(
      variables,
      out
    )
  )
}
