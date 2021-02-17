#' Matrix of Symmetric Paths \eqn{\mathbf{S}}
#'
#' Derives the matrix of symmetric paths (double-headed arrows)
#' \eqn{\mathbf{S}} using the Reticular Action Model (RAM) notation.
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
#' @return \eqn{\mathbf{S} = \mathbf{E}^{-1} \mathbf{C}
#'   \left( \mathbf{E}^{-1} \right)^{\mathsf{T}}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#' @param C `t by t` numeric matrix \eqn{\mathbf{C}}.
#'   Model-implied variance-covariance matrix.
#' @export
S <- function(A,
              C,
              ...) {
  UseMethod("S")
}

#' @rdname S
#' @inheritParams IminusA
#' @inheritParams S
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
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
#' @export
S.default <- function(A,
                      C,
                      ...) {
  stopifnot(is.numeric(C))
  stopifnot(
    matrixR::IsSymmetric(
      round(C, digits = 4)
    )
  )
  stopifnot(identical(dim(A), dim(C)))
  IminusA <- IminusA(A)
  CIminusAt <- tcrossprod(
    x = C,
    y = IminusA
  )
  return(
    IminusA %*% CIminusAt
  )
}

#' @rdname S
#' @inheritParams IminusA
#' @inheritParams S
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
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
#' S(Ryacas::ysym(A), C, R = FALSE, format = "ysym")
#' S(Ryacas::ysym(A), C, R = FALSE, format = "str")
#' S(Ryacas::ysym(A), C, R = FALSE, format = "tex")
#' S(Ryacas::ysym(A), C, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' S(Ryacas::ysym(A), C, R = FALSE, format = "ysym")
#' S(Ryacas::ysym(A), C, R = FALSE, format = "str")
#' S(Ryacas::ysym(A), C, R = FALSE, format = "tex")
#' S(Ryacas::ysym(A), C, R = TRUE)
#' eval(S(Ryacas::ysym(A), C, R = TRUE))
#' @export
S.yac_symbol <- function(A,
                         C,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  IminusAysym <- IminusA(
    A = A,
    exe = FALSE
  )
  # IsSymmetric unsafe for nonnumeric input
  Cysym <- matrixR::MatrixCheck(
    A = yacR::as.ysym(C),
    IsSquareMatrix = TRUE
  )
  ADimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        A,
        ")"
      )
    )
  )
  CDimensions <- as.numeric(
    Ryacas::yac_str(
      paste0(
        "Length(",
        Cysym,
        ")"
      )
    )
  )
  stopifnot(
    identical(
      ADimensions,
      CDimensions
    )
  )
  expr <- paste0(
    IminusAysym,
    "*",
    Cysym,
    "*",
    "Transpose(",
    IminusAysym,
    ")"
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
