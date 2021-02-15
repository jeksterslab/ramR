#' Matrix of Total Effects \eqn{\mathbf{E}}
#'
#' Derives the matrix of total effects \eqn{\mathbf{E}}.
#'
#' The matrix of total effects \eqn{\mathbf{E}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{E}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @return \eqn{\mathbf{E} = \left( \mathbf{I} - \mathbf{A} \right)^{-1}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#'
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' E(A)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' E(Ryacas::ysym(A))
#' E(Ryacas::ysym(A), tex = TRUE)
#' E(Ryacas::ysym(A), ysym = FALSE)
#' E(Ryacas::ysym(A), str = FALSE)
#'
#' beta <- 1
#' E(Ryacas::ysym(A))
#' E(Ryacas::ysym(A), tex = TRUE)
#' E(Ryacas::ysym(A), ysym = FALSE)
#' E(Ryacas::ysym(A), str = FALSE)
#' eval(E(Ryacas::ysym(A), str = FALSE))
#' @export
E <- function(A,
              ...) {
  UseMethod("E")
}

#' @rdname E
#' @inheritParams IminusA
#' @export
E.default <- function(A,
                      ...) {
  return(
    solve(
      IminusA(A)
    )
  )
}

#' @rdname E
#' @inheritParams IminusA
#' @export
E.yac_symbol <- function(A,
                         exe = TRUE,
                         str = TRUE,
                         ysym = TRUE,
                         simplify = FALSE,
                         tex = FALSE,
                         ...) {
  IminusAsym <- IminusA(
    A = A,
    exe = FALSE
  )
  expr <- paste0(
    "Inverse(",
    IminusAsym,
    ")"
  )
  if (exe) {
    return(
      YacExe(
        expr = expr,
        str = str,
        ysym = ysym,
        tex = tex,
        simplify = simplify
      )
    )
  } else {
    return(expr)
  }
}
