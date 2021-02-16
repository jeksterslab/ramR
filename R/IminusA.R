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
#' @export
IminusA <- function(A,
                    ...) {
  UseMethod("IminusA")
}

#' @rdname IminusA
#' @inheritParams IminusA
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' IminusA(A)
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
#' @inheritParams yacR::Exe
#' @param exe Logical.
#'   If `exe = TRUE`,
#'   executes the resulting `yacas` expression.
#'   If `exe = FALSE`,
#'   returns the resulting `yacas` expression as a character string.
#'   If `exe = FALSE`,
#'   the arguments `str`, `ysym`, `simplify`, and `tex`, are ignored.
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "ysym")
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "str")
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "tex")
#' IminusA(Ryacas::ysym(A), R = TRUE)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "ysym")
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "str")
#' IminusA(Ryacas::ysym(A), R = FALSE, format = "tex")
#' IminusA(Ryacas::ysym(A), R = TRUE)
#' eval(IminusA(Ryacas::ysym(A), R = TRUE))
#' @export
IminusA.yac_symbol <- function(A,
                               exe = TRUE,
                               R = FALSE,
                               format = "ysym",
                               simplify = FALSE,
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
      yacR::Exe(
        expr,
        R = R,
        format = format,
        simplify = simplify
      )
    )
  } else {
    return(expr)
  }
}
