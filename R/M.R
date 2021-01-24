#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Covariance Expectations of Observed Variables
#'   \eqn{\mathbf{M}} - Numeric
#'
#' @description Derives the matrix of covariance expectations of observed variables
#'   \eqn{\mathbf{M}}.
#'
#'   \deqn{
#'     \mathbf{M}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{E}
#'     \mathbf{S}
#'     \mathbf{E}^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{C}
#'     \mathbf{F}^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{p \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{p} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{p + q} .
#'
#' @keywords M
#' @family M functions
#' @inheritParams C_num
#' @inherit ramR references
#' @param filter `p x t` numeric matrix
#'   \eqn{\mathbf{F}}.
#'   Filter matrix used to select observed variables.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, 1, 1)
#' diag(S) <- c(0, 0.25, 1)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' colnames(filter) <- c("y", "x", "e")
#' rownames(filter) <- c("y", "x")
#' M_num(A, S, filter)
#' @export
M_num <- function(A,
                  S,
                  filter) {
  # F * (I - A)^{-1} * S * ((I - A)^{-1})^T * F^T
  return(
    filter %*% tcrossprod(
      x = C_num(
        A = A,
        S = S
      ),
      y = filter
    )
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of Covariance Expectations of Observed Variables
#'   \eqn{\mathbf{M}} - Symbolic
#'
#' @keywords M
#' @family M functions
#' @inherit M_num description details references return
#' @inheritParams M_num
#' @inheritParams E_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigma2x", "sigma2e")
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' M_sym(A, S, filter)
#' @export
M_sym <- function(A,
                  S,
                  filter,
                  simplify = FALSE) {
  filter <- Ryacas::ysym(filter)
  out <- filter * C_sym(A, S) * t(filter)
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(
    out
  )
}
