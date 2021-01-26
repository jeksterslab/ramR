#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Total Effects
#'   \eqn{\mathbf{E}} - Numeric
#'
#' @description Derives the matrix of total effects
#'   \eqn{\mathbf{E}}.
#'
#' @details The matrix of total effects
#'   \eqn{\mathbf{E}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{E}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @keywords E
#' @family E functions
#' @inherit ramR references
#' @inheritParams IminusA_num
#' @return Returns the matrix of total effects
#'   \eqn{\mathbf{E}}.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' E_num(A)
#' @export
E_num <- function(A) {
  return(
    solve(IminusA_num(A))
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Total Effects
#'   \eqn{\mathbf{E}} - Symbolic
#'
#' @keywords E
#' @family E functions
#' @inherit E_num description details references return
#' @inheritParams IminusA_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, "beta", 1)
#' E_sym(A)
#' @export
E_sym <- function(A,
                  simplify = FALSE) {
  out <- solve(IminusA_sym(A, simplify))
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(
    out
  )
}
