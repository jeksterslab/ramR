#' Matrix of Covariance Expectations of Observed Variables
#' \eqn{\mathbf{M}}
#'
#' Derives the matrix of covariance expectations
#' of observed variables \eqn{\mathbf{M}}.
#'
#' The matrix of covariance expectations
#' for given variables \eqn{\mathbf{M}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{M}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{E}
#'     \mathbf{S}
#'     \mathbf{E}^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{C}
#'     \mathbf{F}^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{p \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{p} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables,
#'     that is \eqn{p + q} .
#'
#' @return \eqn{
#'   \mathbf{M} = \mathbf{F} \mathbf{C} \mathbf{F}^{\mathsf{T}}
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams C
#' @inheritParams CheckRAMMatrices
#' @export
M <- function(A,
              S,
              Filter = NULL,
              check = TRUE,
              ...) {
  UseMethod("M")
}

#' @rdname M
#' @inheritParams IminusA
#' @inheritParams M
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
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(Filter) <- c("y", "x", "e")
#' M(A, S, Filter)
#' @export
M.default <- function(A,
                      S,
                      Filter = NULL,
                      check = TRUE,
                      ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S,
      Filter = Filter
    )
    A <- RAM$A
    S <- RAM$S
    Filter <- RAM$Filter
  }
  C <- C(
    A = A,
    S = S,
    check = FALSE
  )
  if (is.null(Filter)) {
    return(C)
  } else {
    return(
      Filter %*% tcrossprod(
        x = C,
        y = Filter
      )
    )
  }
}

#' @rdname M
#' @inheritParams IminusA
#' @inheritParams M
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
#' M(Ryacas::ysym(A), S, Filter)
#' M(Ryacas::ysym(A), S, Filter, format = "str")
#' M(Ryacas::ysym(A), S, Filter, format = "tex")
#' M(Ryacas::ysym(A), S, Filter, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#'
#' M(Ryacas::ysym(A), S, Filter)
#' M(Ryacas::ysym(A), S, Filter, format = "str")
#' M(Ryacas::ysym(A), S, Filter, format = "tex")
#' M(Ryacas::ysym(A), S, Filter, R = TRUE)
#' eval(M(Ryacas::ysym(A), S, Filter, R = TRUE))
#' @export
M.yac_symbol <- function(A,
                         S,
                         Filter = NULL,
                         check = TRUE,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  if (check) {
    RAM <- CheckRAMMatrices(A = A, S = S, Filter = Filter)
    S <- RAM$S
    Filter <- RAM$Filter
  } else {
    S <- yacR::as.ysym.mat(S)
    if (!is.null(Filter)) {
      Filter <- yacR::as.ysym.mat(Filter)
    }
  }
  C <- C(
    A = A,
    S = S,
    check = FALSE,
    exe = FALSE
  )
  if (is.null(Filter)) {
    expr <- C
  } else {
    expr <- paste0(
      Filter,
      "*",
      C,
      "*",
      "Transpose(",
      Filter,
      ")"
    )
  }
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
