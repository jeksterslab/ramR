
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Symmetric Paths
#'   \eqn{\mathbf{S}} - Numeric
#'
#' @description Derives the matrix of symmetric paths (double-headed arrows)
#'   \eqn{\mathbf{S}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The matrix of symmetric paths (double-headed arrows)
#'   \eqn{\mathbf{S}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{S}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)
#'     \mathbf{C}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}} \\
#'     =
#'     \mathbf{E}^{-1}
#'     \mathbf{C}
#'     \left(
#'       \mathbf{E}^{-1}
#'     \right)^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{C}_{t \times t}} represents
#'     the model-implied variance-covariance matrix, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @keywords S
#' @family S functions
#' @inheritParams E_num
#' @inherit ramR references
#' @param C `t x t` numeric matrix \eqn{\mathbf{C}}.
#'   Model-implied variance-covariance matrix.
#' @return Returns the matrix of symmetric paths (double-headed arrows)
#'   \eqn{\mathbf{S}}
#'   derived from the `A` and `C` matrices.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, 1, 1)
#' C <- matrix(
#'   data = c(
#'     1.25, 0.25, 1.00,
#'     0.25, 0.25, 0.00,
#'     1.00, 0.00, 1.00
#'   ),
#'   nrow = dim(A)[1]
#' )
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' S_num(A, C)
#' @export
S_num <- function(A,
                  C) {
  if (!matrixR::is_sqr(A)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  if (!identical(dim(A), dim(C))) {
    stop(
      "`A` and `C` should have the same dimensions."
    )
  }
  # I - A
  IminusA <- diag(nrow(A)) - A
  # C * (I - A)^{T}
  CIminusAt <- tcrossprod(
    x = C,
    y = IminusA
  )
  # (I - A) * C * (I - A)^{T}
  return(
    IminusA %*% CIminusAt
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Symmetric Paths
#'   \eqn{\mathbf{S}} - Symmetric
#'
#' @keywords S
#' @family S functions
#' @inherit S_num description details references return
#' @inheritParams S_num
#' @inheritParams E_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, "beta", 1)
#' C <- matrix(
#'   data = c(
#'     "sigma2x*beta^2+sigma2e", "sigma2x*beta", "sigma2e",
#'     "beta*sigma2x", "sigma2x", 0,
#'     "sigma2e", 0, "sigma2e"
#'   ),
#'   nrow = dim(A)[1]
#' )
#' S_sym(A, C)
#' @export
S_sym <- function(A,
                  C,
                  simplify = FALSE) {
  if (!matrixR::is_sqr(A, chk.num = FALSE)) {
    stop(
      "`A` should be a square matrix."
    )
  }
  if (!identical(dim(A), dim(C))) {
    stop(
      "`A` and `C` should have the same dimensions."
    )
  }
  IminusA <- Ryacas::ysym(diag(dim(A)[1])) - Ryacas::ysym(A)
  out <- IminusA * Ryacas::ysym(C) * t(IminusA)
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(
    out
  )
}
