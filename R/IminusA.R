#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title \eqn{\mathbf{I} - \mathbf{A}} - Numeric
#'
#' @description Derives \eqn{\mathbf{I} - \mathbf{A}}
#'
#' @details
#'   \deqn{
#'     \mathbf{I} - \mathbf{A}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @keywords IminusA
#' @family IminusA functions
#' @inherit ramR references
#' @param A `t x t` matrix \eqn{\mathbf{A}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @return Returns \eqn{\mathbf{I} - \mathbf{A}}.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' IminusA_num(A)
#' @export
IminusA_num <- function(A) {
  if (!matrixR::is_sqr(A)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  if (!matrixR::is_nilpot(A)) {
    stop(
      "`A` should be a nilpotent matrix."
    )
  }
  out <- matrixR::ones_from(A) - A
  if (matrixR::is_sing(out)) {
    stop(
      "`I - A` is singular."
    )
  }
  return(out)
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title \eqn{\mathbf{I} - \mathbf{A}} - Symbolic
#'
#' @keywords IminusA
#' @family IminusA functions
#' @inherit IminusA_num description details references return
#' @inheritParams IminusA_num
#' @param simplify Logical.
#'   Simplify the results.
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, "beta", 1)
#' IminusA_sym(A)
#' @export
IminusA_sym <- function(A,
                        simplify = FALSE) {
  if (!matrixR::is_sqr(A, chk.num = FALSE)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  out <- Ryacas::ysym(matrixR::ones_from(A)) - Ryacas::ysym(A)
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(out)
}
