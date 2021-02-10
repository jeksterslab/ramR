#' Matrix of Symmetric Paths \eqn{\mathbf{S}}
#'
#' Derives the matrix of symmetric paths (double-headed arrows) \eqn{\mathbf{S}}
#' using the Reticular Action Model (RAM) notation.
#'
#' The matrix of symmetric paths (double-headed arrows) \eqn{\mathbf{S}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{S}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)
#'     \mathbf{C}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}} \\
#'     =
#'     \mathbf{E}^{-1}
#'     \mathbf{C}
#'     \left(
#'       \mathbf{E}^{-1}
#'     \right)^{\mathsf{T}}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{C}_{t \times t}} represents
#'     the model-implied variance-covariance matrix, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inheritParams IminusA
#' @param C `t x t` numeric matrix \eqn{\mathbf{C}}.
#'   Model-implied variance-covariance matrix.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' C <- matrix(
#'   data = c(
#'     1.25, 0.25, 1.00,
#'     0.25, 0.25, 0.00,
#'     1.00, 0.00, 1.00
#'   ),
#'   nrow = dim(A)[1]
#' )
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' S(A, C)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' C <- matrix(
#'   data = c(
#'     "sigmax2*beta^2+sigmae2", "sigmax2*beta", "sigmae2",
#'     "sigmax2*beta", "sigmax2", 0,
#'     "sigmae2", 0, "sigmae2"
#'   ),
#'   nrow = dim(A)[1]
#' )
#' S(Ryacas::ysym(A), C)
#' S(Ryacas::ysym(A), C, tex = TRUE)
#' S(Ryacas::ysym(A), C, ysym = FALSE)
#' S(Ryacas::ysym(A), C, str = FALSE)
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' S(Ryacas::ysym(A), C)
#' S(Ryacas::ysym(A), C, tex = TRUE)
#' S(Ryacas::ysym(A), C, ysym = FALSE)
#' S(Ryacas::ysym(A), C, str = FALSE)
#' eval(S(Ryacas::ysym(A), C, str = FALSE))
#' @export
S <- function(A,
              C,
              ...) {
  UseMethod("S")
}

#' @rdname S
#' @inheritParams IminusA
#' @inheritParams S
#' @export
S.default <- function(A,
                      C,
                      ...) {
  stopifnot(matrixR::IsSymmetric(C))
  stopifnot(identical(dim(A), dim(C)))
  # I - A
  IminusA <- IminusA.default(A)
  # C * (I - A)^{T}
  CIminusAt <- tcrossprod(
    x = C,
    y = IminusA
  )
  # (I - A) * C * (I - A)^{T}
  return(
    IminusA %*% CIminusAt
  )
}

#' @rdname S
#' @inheritParams IminusA
#' @inheritParams S
#' @export
S.yac_symbol <- function(A,
                         C,
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
  if (isFALSE(methods::is(C, "yac_symbol"))) {
    C <- Ryacas::ysym(C)
  }
  x_res <- Ryacas::yac_str(C$yacas_cmd)
  x <- Ryacas::ysym(x_res)
  stopifnot(x$is_mat)
  stopifnot(matrixR::IsSymmetric(x))
  At <- as.numeric(Ryacas::yac_str(paste0("Length(", y, ")")))
  Ct <- as.numeric(Ryacas::yac_str(paste0("Length(", x, ")")))
  if (isFALSE(identical(At, Ct))) {
    stop(
      "`A` and `C` do not have the same dimensions."
    )
  }
  I <- paste0("Identity(Length(", y, "))")
  expr <- paste0(
    "(", I, "-", y, ")", "*", x, "*",
    "Transpose(", I, "-", y, ")"
  )
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
