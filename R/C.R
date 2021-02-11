#' Matrix of Covariance Expectations \eqn{\mathbf{C}}
#'
#' Derives the matrix of covariance expectations \eqn{\mathbf{C}}.
#'
#' The matrix of covariance expectations \eqn{\mathbf{C}}
#' as a function of Reticular Action Model (RAM) matrices is given by
#'
#'   \deqn{
#'     \mathbf{C}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\mathsf{T}} \\
#'     =
#'     \mathbf{E}
#'     \mathbf{S}
#'     \mathbf{E}^{\mathsf{T}}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances, and
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inheritParams IminusA
#' @param S `t x t` numeric matrix \eqn{\mathbf{S}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' diag(S) <- c(0, 0.25, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' C(A, S)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigmax2", "sigmae2")
#' C(Ryacas::ysym(A), S)
#' C(Ryacas::ysym(A), S, tex = TRUE)
#' C(Ryacas::ysym(A), S, ysym = FALSE)
#' C(Ryacas::ysym(A), S, str = FALSE)
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' C(Ryacas::ysym(A), S)
#' C(Ryacas::ysym(A), S, tex = TRUE)
#' C(Ryacas::ysym(A), S, ysym = FALSE)
#' C(Ryacas::ysym(A), S, str = FALSE)
#' eval(C(Ryacas::ysym(A), S, str = FALSE))
#' @export
C <- function(A,
              S,
              ...) {
  UseMethod("C")
}

#' @rdname C
#' @inheritParams IminusA
#' @inheritParams C
#' @export
C.default <- function(A,
                      S,
                      ...) {
  stopifnot(
    matrixR::IsSymmetric(
      S
    )
  )
  stopifnot(
    identical(
      dim(A),
      dim(S)
    )
  )
  # (I - A)^{-1}
  E <- E.default(A)
  # S * ((I - A)^{-1})^T
  STinvIminusA <- tcrossprod(
    x = S,
    y = E
  )
  # (I - A)^{-1} * S * ((I - A)^{-1})^T
  return(
    E %*% STinvIminusA
  )
}

#' @rdname C
#' @inheritParams IminusA
#' @inheritParams C
#' @export
C.yac_symbol <- function(A,
                         S,
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
  # apply IsNilpotent in the future
  if (methods::is(S, "yac_symbol")) {
    Sysym <- S
  } else {
    Sysym <- Ryacas::ysym(
      S
    )
  }
  Sysym <- Ryacas::ysym(
    Ryacas::yac_str(
      Sysym$yacas_cmd
    )
  )
  stopifnot(
    Sysym$is_mat
  )
  stopifnot(
    matrixR::IsSymmetric(
      Sysym
    )
  )
  ADimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        Aysym,
        ")"
      )
    )
  )
  SDimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        Sysym,
        ")"
      )
    )
  )
  stopifnot(
    identical(
      ADimensions,
      SDimensions
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
    Sysym,
    "*",
    "Transpose(",
    E,
    ")"
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
