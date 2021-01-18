#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Expected Values
#'   \eqn{\mathbf{v}} - Numeric
#'
#' @description Derives the vector of expected values
#'   \eqn{\mathbf{v}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The vector of expected values
#'   \eqn{\mathbf{v}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{v}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{E}
#'     \mathbf{u}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix, and
#'   - \eqn{\mathbf{u}_{t \times 1}} vector of parameters for the mean structure.
#'
#' @keywords v
#' @family v functions
#' @inheritParams C_num
#' @param u vector of length `t` or `t by 1` matrix.
#'   Mean structure parameters.
#' @return Returns the vector of expected values
#'   \eqn{\mathbf{v}}.
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
#' A[1, 2] <- 1
#' A[1, 3] <- 1
#' u <- c(0.00, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' v_num(A, u)
#' @export
v_num <- function(A,
                  u) {
  out <- as.matrix(
    E_num(
      A
    ) %*% as.matrix(u)
  )
  colnames(out) <- "v"
  return(
    out
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Expected Values
#'   \eqn{\mathbf{v}} - Symbolic
#'
#' @keywords v
#' @family v functions
#' @inherit v_num description details references return
#' @inheritParams v_num
#' @examples
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
#' u <- c("alpha", "mux", 0)
#' v_sym(A, u)
#' @export
v_sym <- function(A,
                  u) {
  u <- as.matrix(u)
  return(
    E_sym(A) * Ryacas::ysym(u)
  )
}
