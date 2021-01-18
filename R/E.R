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
#' @param A `t x t` matrix \eqn{\mathbf{A}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @return Returns the matrix of total effects
#'   \eqn{\mathbf{E}}.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = -0.5 + 1 * x + e
#'
#' A <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, 2] <- 1
#' A[1, 3] <- 1
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' E_num(A)
#' @export
E_num <- function(A) {
  if (!matrixR::is_sqr(A)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  IminusA <- diag(nrow(A)) - A
  if (matrixR::is_sing(IminusA)) {
    stop(
      "`I - A` is singular."
    )
  }
  return(
    solve(IminusA)
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
#' @inheritParams E_num
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, 2] <- "beta"
#' A[1, 3] <- 1
#' E_sym(A)
#' @export
E_sym <- function(A) {
  if (!matrixR::is_sqr(A, chk.num = FALSE)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  return(
    solve(Ryacas::ysym(diag(dim(A)[1])) - Ryacas::ysym(A))
  )
}
