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
#' @return \eqn{
#'   \mathbf{C} = \mathbf{E} \mathbf{S} \mathbf{E}^{\mathsf{T}}
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams IminusA
#' @inheritParams CheckRAMMatrices
#' @export
C <- function(A,
              S,
              check = TRUE,
              ...) {
  UseMethod("C")
}

#' @rdname C
#' @inheritParams IminusA
#' @inheritParams C
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#' #--------------------------------------------------------------------
#'
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' diag(S) <- c(0, 0.25, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' C(A, S)
#' @export
C.default <- function(A,
                      S,
                      check = TRUE,
                      ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S
    )
    A <- RAM$A
    S <- RAM$S
  }
  E <- E(
    A = A,
    check = FALSE
  )
  STinvIminusA <- tcrossprod(
    x = S,
    y = E
  )
  return(
    E %*% STinvIminusA
  )
}

#' @rdname C
#' @inheritParams IminusA
#' @inheritParams C
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#' #--------------------------------------------------------------------
#'
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigmax2", "sigmae2")
#' C(Ryacas::ysym(A), S)
#' C(Ryacas::ysym(A), S, format = "str")
#' C(Ryacas::ysym(A), S, format = "tex")
#' C(Ryacas::ysym(A), S, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#'
#' C(Ryacas::ysym(A), S)
#' C(Ryacas::ysym(A), S, format = "str")
#' C(Ryacas::ysym(A), S, format = "tex")
#' C(Ryacas::ysym(A), S, R = TRUE)
#' eval(C(Ryacas::ysym(A), S, R = TRUE))
#' @export
C.yac_symbol <- function(A,
                         S,
                         check = TRUE,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  if (check) {
    S <- CheckRAMMatrices(
      A = A,
      S = S
    )$S
  } else {
    S <- yacR::as.ysym.mat(S)
  }
  E <- E(
    A = A,
    check = FALSE,
    exe = FALSE
  )
  expr <- paste0(
    E,
    "*",
    S,
    "*",
    "Transpose(",
    E,
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
