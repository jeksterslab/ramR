#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Mean Structure Parameters
#'   \eqn{\mathbf{u}} - Numeric
#'
#' @description Derives the vector of mean structure parameters
#'   \eqn{\mathbf{u}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The vector of mean structure parameters
#'   \eqn{\mathbf{v}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{u}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)
#'     \mathbf{v} \\
#'     =
#'     \mathbf{E}^{-1}
#'     \mathbf{v}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix, and
#'   - \eqn{\mathbf{v}_{t \times 1}} vector of expected values.
#'
#' @keywords u
#' @family u functions
#' @inheritParams C_num
#' @param v vector of length `t` or `t by 1` matrix.
#'   Expected values.
#' @return Returns the vector of mean structire parameters
#'   \eqn{\mathbf{u}}.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, 1, 1)
#' v <- c(0.50, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' u_num(A, v)
#' @export
u_num <- function(A,
                  v) {
  # (I - A) * v
  out <- as.matrix(
    (matrixR::ones_from(A) - A) %*% as.matrix(v)
  )
  colnames(out) <- "u"
  return(
    out
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Mean Structure Parameters
#'   \eqn{\mathbf{u}} - Symbolic
#'
#' @keywords u
#' @family u functions
#' @inherit u_num description details references return
#' @inheritParams u_num
#' @inheritParams IminusA_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, "beta", 1)
#' v <- c("alpha+beta*mux", "mux", 0)
#' u_sym(A, v)
#' @export
u_sym <- function(A,
                  v,
                  simplify = FALSE) {
  v <- as.matrix(v)
  out <- (Ryacas::ysym(matrixR::ones_from(A)) - Ryacas::ysym(A)) * Ryacas::ysym(v)
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(
    out
  )
}
