#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Covariance Expectations
#'   \eqn{\mathbf{C}} - Numeric
#'
#' @description Derives the matrix of covariance expectations
#'   \eqn{\mathbf{C}}.
#'
#' @details The matrix of covariance expectations
#'   \eqn{\mathbf{C}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{C}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\mathsf{T}} \\
#'     =
#'     \mathbf{E}
#'     \mathbf{S}
#'     \mathbf{E}^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @keywords C
#' @family C functions
#' @inheritParams S_num
#' @inherit ramR references
#' @param S `t x t` numeric matrix
#'   \eqn{\mathbf{S}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = -0.5 + 1 * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, 2] <- 1
#' A[1, 3] <- 1
#' diag(S) <- c(0, 0.25, 0.25)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' C_num(A, S)
#' @export
C_num <- function(A,
                  S) {
  if (!matrixR::is_sqr(S)) {
    stop(
      "`S` should be a square matrix."
    )
  }
  if (!identical(dim(A), dim(S))) {
    stop(
      "`A` and `S` should have the same dimensions."
    )
  }
  E <- E_num(
    A = A
  )
  # S * ((I - A)^{-1})^T
  STinvIminusA <- tcrossprod(
    x = S,
    y = E
  )
  # (I - A)^{-1} * S * ((I - A)^{-1})^T
  return(
    E %*% STinvIminusA
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Covariance Expectations
#'   \eqn{\mathbf{C}} - Symbolic
#'
#' @keywords C
#' @family C functions
#' @inherit C_num description details references return
#' @inheritParams C_num
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, 2] <- "beta"
#' A[1, 3] <- 1
#' diag(S) <- c(0, "sigma2x", "sigma2e")
#' C_sym(A, S)
#' @export
C_sym <- function(A,
                  S) {
  if (!matrixR::is_sqr(S, chk.num = FALSE)) {
    stop(
      "`S` should be a square matrix."
    )
  }
  if (!identical(dim(A), dim(S))) {
    stop(
      "`A` and `S` should have the same dimensions."
    )
  }
  E <- E_sym(A)
  return(
    E * Ryacas::ysym(S) * t(E)
  )
}
