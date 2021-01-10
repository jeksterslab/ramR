#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Total Effects \eqn{\mathbf{E}}
#'
#' @description Derives the matrix of total effects \eqn{\mathbf{E}}.
#'
#' @details The matrix of total effects \eqn{\mathbf{E}} is given by
#'   \deqn{
#'     \mathbf{E}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   }
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A `t x t` numeric matrix
#'   \eqn{\mathbf{A}_{t \times t}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @export
E <- function(A) {
  if (!matrixR::is_sqr(X = A)) {
    stop(
      "`A` should be a square matrix."
    )
  } else {
    IminusA <- diag(nrow(A)) - A
    if (matrixR::is_inv(IminusA)) {
      return(solve(IminusA))
    } else {
      stop(
        "`I - A` is not invertible."
      )
    }
  }
}
