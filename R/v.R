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
#'
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' u <- c(0.00, 0.50, 0.00)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' v(A, u)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' u <- c("alpha", "mux", 0)
#' v(Ryacas::ysym(A), u)
#' v(Ryacas::ysym(A), u, tex = TRUE)
#' v(Ryacas::ysym(A), u, ysym = FALSE)
#' v(Ryacas::ysym(A), u, str = FALSE)
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' v(Ryacas::ysym(A), u)
#' v(Ryacas::ysym(A), u, tex = TRUE)
#' v(Ryacas::ysym(A), u, ysym = FALSE)
#' v(Ryacas::ysym(A), u, str = FALSE)
#' eval(v(Ryacas::ysym(A), u, str = FALSE))
#' @export
v <- function(A,
              u,
              ...) {
  UseMethod("v")
}

#' @rdname v
#' @inheritParams IminusA
#' @inheritParams v
#' @export
v.default <- function(A,
                      u,
                      ...) {
  u <- matrix(
    u,
    ncol = 1
  )
  stopifnot(identical(dim(A)[1], dim(u)[1]))
  E <- E.default(A)
  v <- as.matrix(E %*% u)
  colnames(v) <- "v"
  return(v)
}

#' @rdname v
#' @inheritParams IminusA
#' @inheritParams v
#' @export
v.yac_symbol <- function(A,
                         u,
                         str = TRUE,
                         ysym = TRUE,
                         simplify = FALSE,
                         tex = FALSE,
                         ...) {
  stopifnot(
    methods::is(
      A,
      "yac_symbol"
    )
  )
  Aysym <- Ryacas::ysym(
    Ryacas::yac_str(
      A$yacas_cmd
    )
  )
  stopifnot(
    Aysym$is_mat
  )
  stopifnot(
    matrixR::IsSquareMatrix(
      Aysym
    )
  )
  I <- paste0(
    "Identity(Length(",
    Aysym,
    "))"
  )
  E <- paste0(
    "Inverse(",
    I,
    "-",
    Aysym,
    ")"
  )
  expr <- paste0(
    E,
    "*",
    Ryacas::ysym(
      matrix(
        u,
        ncol = 1
      )
    )
  )
  return(
    YacExe(
      expr = expr,
      str = str,
      ysym = ysym,
      simplify = simplify,
      tex = tex
    )
  )
}
