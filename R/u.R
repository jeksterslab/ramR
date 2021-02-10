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
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inheritParams IminusA
#' @param v vector of length `t` or `t by 1` matrix.
#'   Expected values.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' v <- c(0.50, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' u(A, v)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' v <- c("alpha+beta*mux", "mux", 0)
#' # u(A, v, sym = TRUE)
#' # u(A, v, sym = TRUE, str = FALSE)
#' # u(A, v, sym = TRUE, ysym = FALSE)
#' # u(A, v, sym = TRUE, tex = TRUE)
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' u(Ryacas::ysym(A), v)
#' u(Ryacas::ysym(A), v, tex = TRUE)
#' u(Ryacas::ysym(A), v, ysym = FALSE)
#' u(Ryacas::ysym(A), v, str = FALSE)
#' eval(u(Ryacas::ysym(A), v, str = FALSE))
#' @export
u <- function(A,
              v,
              ...) {
  UseMethod("u")
}

#' @rdname u
#' @inheritParams IminusA.default
#' @inheritParams u
#' @export
u.default <- function(A,
                      v,
                      ...) {
  v <- matrix(
    v,
    ncol = 1
  )
  stopifnot(identical(dim(A)[1], dim(v)[1]))
  IminusA <- IminusA.default(A)
  u <- as.matrix(IminusA %*% v)
  colnames(u) <- "u"
  return(u)
}

#' @rdname u
#' @inheritParams IminusA.yac_symbol
#' @inheritParams u
#' @export
u.yac_symbol <- function(A,
                         v,
                         str = TRUE,
                         ysym = TRUE,
                         simplify = FALSE,
                         tex = FALSE,
                         ...) {
  stopifnot(methods::is(A, "yac_symbol"))
  y_res <- Ryacas::yac_str(A$yacas_cmd)
  y <- Ryacas::ysym(y_res)
  stopifnot(y$is_mat)
  stopifnot(matrixR::IsSquareMatrix(y))
  I <- paste0("Identity(Length(", y, "))")
  v <- matrix(
    v,
    ncol = 1
  )
  a <- Ryacas::ysym(v)
  expr <- paste0(I, "*", a)
  return(
    .exe(
      expr = expr,
      str = str,
      ysym = ysym,
      simplify = simplify,
      tex = tex
    )
  )
}
