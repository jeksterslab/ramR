#' Equations to Expectations
#'
#' Converts equations to expectations.
#' The argument `eq` is a character string
#' that specifies the associations between the variables.
#' See `Syntax`, `Operations`, `par.label`,
#' `Line breaks`, and `Comments` below.
#'
#' @inherit Expectations details
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs <space> op <space> rhs <space> par.label <\n> or <;>`
#'
#'   \describe{
#'     \item{lhs}{
#'       is the variable on the **left-hand side**,
#'     }
#'     \item{rhs}{
#'       is the variable on the **right-hand side**,
#'     }
#'     \item{op}{
#'       is the **operation** between `lhs` and `rhs`,
#'     }
#'     \item{par.label}{
#'       is the column of **parameter label**,
#'     }
#'     \item{\\n or ;}{
#'       are **line breaks**.
#'       **Each line should end with a line break.**
#'     }
#'   }
#'
#' @section Operations:
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{
#'       `left-hand side` measured **by** `right-hand side`,
#'     }
#'     \item{on}{
#'       `left-hand side` regressed **on** `right-hand side`,
#'     }
#'     \item{with}{
#'       `left-hand side` covarying **with** `right-hand side`,
#'     }
#'     \item{on 1}{
#'       `left-hand side` regressed **on 1** for mean structure.
#'     }
#'   }
#'
#' @section par.label:
#'   Each parameter should be labeled.
#'   The `par.label` should be a number for fixed parameters
#'   and a character string for free parameters.
#'   Equality contraints can be imposed by using the same `par.label`.
#'
#' @section Line breaks:
#'   The characters `\n` and `;` can be used as line breaks.
#'   **Each line should end with a line break.**
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @return Returns a list with the following elements
#'
#'   \describe{
#'     \item{par.table}{Parameter table.}
#'     \item{variables}{Variable names.}
#'     \item{g.variables}{Variable names of observed variables.}
#'     \item{h.variables}{Variable names of latent variables.}
#'     \item{A}{
#'       `t by t` matrix \eqn{\mathbf{A}}.
#'       Asymmetric paths (single-headed arrows),
#'       such as regression coefficients and factor loadings.
#'     }
#'     \item{S}{
#'       `t by t` numeric matrix \eqn{\mathbf{S}}.
#'       Symmetric paths (double-headed arrows),
#'       such as variances and covariances.
#'     }
#'     \item{u}{
#'       `t by 1` matrix \eqn{\mathbf{u}} of mean structure parameters.
#'     }
#'     \item{Filter}{
#'       `p by t` numeric matrix
#'       \eqn{\mathbf{F}}.
#'       Filter matrix used to select observed variables.
#'     }
#'     \item{v}{
#'       `t by 1` matrix \eqn{\mathbf{v}}
#'       of expected values.
#'     }
#'     \item{g}{
#'       `p by 1` matrix \eqn{\mathbf{g}}
#'       of expected values of observed variables.
#'     }
#'     \item{C}{
#'       `t by t` matrix \eqn{\mathbf{C}}
#'       of expected covariances.
#'     }
#'     \item{M}{
#'       `p by p` matrix \eqn{\mathbf{M}}
#'       of expected covariances of observed variables.
#'     }
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family eq functions
#' @keywords eq
#'
#' @inherit ramR references
#'
#' @inheritParams Eq2RAM
#' @inheritParams Expectations
#'
#' @examples
#' # Numerical ---------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs par.label
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' Eq2Expectations(eq)
#'
#' # Symbolic ----------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs par.label
#'     e   by   y   1
#'     y   on   x   beta
#'     e   with e   sigmae2
#'     x   with x   sigmax2
#'     y   on   1   alpha
#'     x   on   1   mux
#' "
#' Eq2Expectations(eq, par = FALSE)
#' Eq2Expectations(eq, par = TRUE)
#'
#' # Expressions using `par.label`
#'
#' beta <- 1
#' sigmae2 <- 1
#' sigmax2 <- 0.25
#' alpha <- 0
#' mux <- 0.50
#' Exp <- Eq2Expectations(eq, par = FALSE, R = TRUE)
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
#' Exp <- Eq2Expectations(eq, par = TRUE, R = TRUE)
#' eval(Exp$M)
#' eval(Exp$g)
#' @export
Eq2Expectations <- function(eq,
                            par = FALSE,
                            check = TRUE,
                            R = FALSE,
                            format = "ysym",
                            simplify = FALSE) {
  Expectations <- Eq2RAM(
    eq,
    par = par
  )
  if (is.numeric(Expectations$par.table$par.label)) {
    out <- Expectations.default(
      A = Expectations$A,
      S = Expectations$S,
      u = Expectations$u,
      Filter = Expectations$Filter,
      check = check
    )
  } else {
    out <- Expectations.yac_symbol(
      A = Ryacas::ysym(Expectations$A),
      S = Expectations$S,
      u = Expectations$u,
      Filter = Expectations$Filter,
      check = check,
      R = R,
      format = format,
      simplify = simplify
    )
  }
  return(
    c(
      par.table = list(Expectations[["par.table"]]),
      variables = list(Expectations[["variables"]]),
      g.variables = list(Expectations[["g.variables"]]),
      h.variables = list(Expectations[["h.variables"]]),
      out
    )
  )
}
