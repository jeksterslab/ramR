#' \eqn{\mathbf{I} - \mathbf{A}}
#'
#' Derives \eqn{\mathbf{I} - \mathbf{A}}
#'
#'   \deqn{
#'     \mathbf{I} - \mathbf{A}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @return \eqn{\mathbf{I} - \mathbf{A}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @param A `t by t` matrix \eqn{\mathbf{A}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param ... ...
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
#' IminusA(A)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' IminusA(Ryacas::ysym(A))
#' IminusA(Ryacas::ysym(A), ysym = FALSE)
#' IminusA(Ryacas::ysym(A), tex = TRUE)
#' IminusA(Ryacas::ysym(A), str = FALSE)
#'
#' beta <- 1
#' IminusA(Ryacas::ysym(A))
#' IminusA(Ryacas::ysym(A), ysym = FALSE)
#' IminusA(Ryacas::ysym(A), tex = TRUE)
#' IminusA(Ryacas::ysym(A), str = FALSE)
#' eval(IminusA(Ryacas::ysym(A), str = FALSE))
#' @export
IminusA <- function(A,
                    ...) {
  UseMethod("IminusA")
}

#' @rdname IminusA
#' @inheritParams IminusA
#' @export
IminusA.default <- function(A,
                            ...) {
  stopifnot(is.numeric(A))
  stopifnot(
    matrixR::IsNilpotent(
      A
    )
  )
  return(matrixR::IdentityFrom(A) - A)
}

#' @rdname IminusA
#' @inheritParams IminusA
#' @inheritParams YacExe
#' @param exe Logical.
#'   If `exe = TRUE`, executes the resulting `yacas` expression.
#'   If `exe = FALSE`, returns the resulting `yacas` expression as a character string.
#'   If `exe = FALSE`, the arguments `str`, `ysym`, `simplify`, and `tex`, are ignored.
#' @export
IminusA.yac_symbol <- function(A,
                               exe = TRUE,
                               str = TRUE,
                               ysym = TRUE,
                               simplify = FALSE,
                               tex = FALSE,
                               ...) {
  Aysym <- matrixR::MatrixCheck(
    A = A,
    IsSquareMatrix = TRUE
  )
  expr <- paste0(
    "Identity(Length(",
    Aysym,
    "))",
    "-",
    Aysym
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
