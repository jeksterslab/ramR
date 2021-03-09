#' Matrix of Total Effects \eqn{\mathbf{E}}
#'
#' Derives the matrix of total effects \eqn{\mathbf{E}}.
#'
#' The matrix of total effects \eqn{\mathbf{E}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{E}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @return \eqn{
#'   \mathbf{E} = \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#' @export
E <- function(A,
              check = TRUE,
              ...) {
  UseMethod("E")
}

#' @rdname E
#' @inheritParams IminusA
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#' #--------------------------------------------------------------------
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' E(A)
#' @export
E.default <- function(A,
                      check = TRUE,
                      ...) {
  return(
    solve(
      IminusA(
        A = A,
        check = check
      )
    )
  )
}

#' @rdname E
#' @inheritParams IminusA
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#' #--------------------------------------------------------------------
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' E(Ryacas::ysym(A))
#' E(Ryacas::ysym(A), format = "str")
#' E(Ryacas::ysym(A), format = "tex")
#' E(Ryacas::ysym(A), R = TRUE)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#'
#' E(Ryacas::ysym(A))
#' E(Ryacas::ysym(A), format = "str")
#' E(Ryacas::ysym(A), format = "tex")
#' E(Ryacas::ysym(A), R = TRUE)
#' eval(E(Ryacas::ysym(A), R = TRUE))
#' @export
E.yac_symbol <- function(A,
                         check = TRUE,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  expr <- paste0(
    "Inverse(",
    IminusA(
      A = A,
      check = check,
      exe = FALSE
    ),
    ")"
  )
  if (exe) {
    return(
      yacR::Exe(
        expr,
        R = R,
        format = format,
        simplify = simplify
      )
    )
  } else {
    return(expr)
  }
}
