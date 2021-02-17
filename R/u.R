#' Vector of Mean Structure Parameters \eqn{\mathbf{u}}
#'
#' Derives the vector of mean structure parameters \eqn{\mathbf{u}}
#' using the Reticular Action Model (RAM) notation.
#'
#' The vector of mean structure parameters \eqn{\mathbf{v}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
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
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix, and
#'   - \eqn{\mathbf{v}_{t \times 1}} vector of expected values.
#'
#' @return \eqn{\mathbf{u} = \mathbf{E}^{-1} \mathbf{v}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#' @param v vector of length `t` or `t by 1` matrix.
#'   Expected values.
#' @export
u <- function(A,
              v,
              ...) {
  UseMethod("u")
}

#' @rdname u
#' @inheritParams IminusA.default
#' @inheritParams u
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' v <- c(0.50, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' u(A, v)
#' @export
u.default <- function(A,
                      v,
                      ...) {
  IminusA <- IminusA(A)
  stopifnot(is.numeric(v))
  v <- as.matrix(v)
  stopifnot(identical(dim(A)[1], dim(v)[1]))
  u <- as.matrix(IminusA %*% v)
  colnames(u) <- "u"
  return(u)
}

#' @rdname u
#' @inheritParams IminusA.yac_symbol
#' @inheritParams u
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' v <- c("alpha+beta*mux", "mux", 0)
#' u(Ryacas::ysym(A), v)
#' u(Ryacas::ysym(A), v, format = "tex")
#' u(Ryacas::ysym(A), v, format = "ysym")
#' u(Ryacas::ysym(A), v, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' u(Ryacas::ysym(A), v)
#' u(Ryacas::ysym(A), v, format = "tex")
#' u(Ryacas::ysym(A), v, format = "ysym")
#' u(Ryacas::ysym(A), v, R = TRUE)
#' eval(u(Ryacas::ysym(A), v, R = TRUE))
#' @export
u.yac_symbol <- function(A,
                         v,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  Aysym <- matrixR::MatrixCheck(
    A,
    IsSquareMatrix = TRUE
  )
  Iysym <- paste0(
    "Identity(Length(",
    Aysym,
    "))"
  )
  vysym <- yacR::as.ysym.mat(v)
  ADimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        Aysym,
        ")"
      )
    )
  )
  vDimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        vysym,
        ")"
      )
    )
  )
  stopifnot(
    identical(
      ADimensions,
      vDimensions
    )
  )
  expr <- paste0(
    Iysym,
    "*",
    vysym
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
