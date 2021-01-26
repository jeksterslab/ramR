#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Expected Values of Observed Variables
#'   \eqn{\mathbf{g}} - Numeric
#'
#' @description Derives the vector of expected values of observed variables
#'   \eqn{\mathbf{g}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The vector of expected values of observed variables
#'   \eqn{\mathbf{g}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \mathbf{g}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{E}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{v}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\mathbf{u}_{t \times 1}} vector of parameters for the mean structurem,
#'   - \eqn{\mathbf{F}_{p \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{p} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{p + q} .
#'
#' @keywords g
#' @family g functions
#' @inheritParams M_num
#' @param u vector of length `t` or `t by 1` matrix.
#'   Mean structure parameters.
#' @return Returns the vector of expected values of observed variables
#'   \eqn{\mathbf{g}}.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, 1, 1)
#' u <- c(0.00, 0.50, 0.00)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' colnames(filter) <- c("y", "x", "e")
#' rownames(filter) <- c("y", "x")
#' g_num(A, u, filter)
#' @export
g_num <- function(A,
                  u,
                  filter) {
  out <- as.matrix(
    filter %*% v_num(
      A,
      u
    )
  )
  colnames(out) <- "g"
  return(
    out
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Vector of Expected Values of Observed Variables
#'   \eqn{\mathbf{g}} - Symbolic
#'
#' @keywords g
#' @family g functions
#' @inherit g_num description details references return
#' @inheritParams g_num
#' @inheritParams IminusA_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- matrixR::zeroes(3, 3)
#' A[1, ] <- c(0, "beta", 1)
#' u <- c("alpha", "mux", 0)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' g_sym(A, u, filter)
#' @export
g_sym <- function(A,
                  u,
                  filter,
                  simplify = FALSE) {
  out <- Ryacas::ysym(filter) * v_sym(
    A,
    u
  )
  if (simplify) {
    out <- Ryacas::simplify(out)
  }
  return(
    out
  )
}
