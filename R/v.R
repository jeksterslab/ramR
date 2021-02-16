#' Vector of Expected Values \eqn{\mathbf{v}}
#'
#' Derives the vector of expected values \eqn{\mathbf{v}}
#' using the Reticular Action Model (RAM) notation.
#'
#' The vector of expected values \eqn{\mathbf{v}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
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
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix, and
#'   - \eqn{\mathbf{u}_{t \times 1}} vector of parameters
#'     for the mean structure.
#'
#' @return \eqn{\mathbf{v} = \mathbf{E} \mathbf{u}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#' @param u vector of length `t` or `t by 1` matrix.
#'   Mean structure parameters.
#' @export
v <- function(A,
              u,
              ...) {
  UseMethod("v")
}

#' @rdname v
#' @inheritParams IminusA
#' @inheritParams v
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' u <- c(0.00, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' v(A, u)
#' @export
v.default <- function(A,
                      u,
                      ...) {
  E <- E(A)
  stopifnot(is.numeric(u))
  u <- as.matrix(u)
  stopifnot(identical(dim(A)[1], dim(u)[1]))
  v <- as.matrix(E %*% u)
  colnames(v) <- "v"
  return(v)
}

#' @rdname v
#' @inheritParams IminusA
#' @inheritParams v
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' u <- c("alpha", "mux", 0)
#' v(Ryacas::ysym(A), u)
#' v(Ryacas::ysym(A), u, format = "tex")
#' v(Ryacas::ysym(A), u, format = "ysym")
#' v(Ryacas::ysym(A), u, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' v(Ryacas::ysym(A), u)
#' v(Ryacas::ysym(A), u, format = "tex")
#' v(Ryacas::ysym(A), u, format = "ysym")
#' v(Ryacas::ysym(A), u, R = TRUE)
#' eval(v(Ryacas::ysym(A), u, R = TRUE))
#' @export
v.yac_symbol <- function(A,
                         u,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  Eysym <- E(
    A = A,
    exe = FALSE
  )
  uysym <- RVector2Yac(u, col = TRUE)
  expr <- paste0(
    Eysym,
    "*",
    uysym
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
