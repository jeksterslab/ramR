#' Matrix of Covariance Expectations of Observed Variables
#' \eqn{\mathbf{M}}
#'
#' Derives the matrix of covariance expectations of observed variables
#' \eqn{\mathbf{M}}.
#'
#' The matrix of covariance expectations for given variables
#' \eqn{\mathbf{M}} as a function of Reticular Action Model (RAM) matrices
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
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inheritParams C
#' @param Filter `p x t` numeric matrix
#'   \eqn{\mathbf{F}}.
#'   Filter matrix used to select observed variables.
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
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(Filter) <- c("y", "x", "e")
#' M(A, S, Filter)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigmax2", "sigmae2")
#' M(Ryacas::ysym(A), S, Filter)
#' M(Ryacas::ysym(A), S, Filter, tex = TRUE)
#' M(Ryacas::ysym(A), S, Filter, ysym = FALSE)
#' M(Ryacas::ysym(A), S, Filter, str = FALSE)
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' M(Ryacas::ysym(A), S, Filter)
#' M(Ryacas::ysym(A), S, Filter, tex = TRUE)
#' M(Ryacas::ysym(A), S, Filter, ysym = FALSE)
#' M(Ryacas::ysym(A), S, Filter, str = FALSE)
#' eval(M(Ryacas::ysym(A), S, Filter, str = FALSE))
#' @export
M <- function(A,
              S,
              Filter,
              ...) {
  UseMethod("M")
}

#' @rdname M
#' @inheritParams IminusA
#' @inheritParams M
#' @export
M.default <- function(A,
                      S,
                      Filter = NULL,
                      ...) {
  C <- C.default(A, S)
  if (isFALSE(is.null(Filter))) {
    if (dim(C)[1] != dim(Filter)[2]) {
      stop(
        "`A` and `Filter` do not have compatible dimensions."
      )
    }
    return(
      Filter %*% tcrossprod(
        x = C,
        y = Filter
      )
    )
  }
  return(C)
}

#' @rdname M
#' @inheritParams IminusA
#' @inheritParams M
#' @export
M.yac_symbol <- function(A,
                         S,
                         Filter = NULL,
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
  if (isFALSE(methods::is(S, "yac_symbol"))) {
    S <- Ryacas::ysym(S)
  }
  x_res <- Ryacas::yac_str(S$yacas_cmd)
  x <- Ryacas::ysym(x_res)
  stopifnot(x$is_mat)
  stopifnot(matrixR::IsSymmetric(x))
  At <- as.numeric(Ryacas::yac_str(paste0("Length(", y, ")")))
  St <- as.numeric(Ryacas::yac_str(paste0("Length(", x, ")")))
  if (isFALSE(identical(At, St))) {
    stop(
      "`A` and `S` do not have the same dimensions."
    )
  }
  I <- paste0("Identity(Length(", y, "))")
  E <- paste0("Inverse(", I, "-", y, ")")
  C <- paste0(E, "*", x, "*", "Transpose(", E, ")")
  if (isFALSE(is.null(Filter))) {
    if (isFALSE(methods::is(Filter, "yac_symbol"))) {
      Filter <- Ryacas::ysym(Filter)
    }
    z_res <- Ryacas::yac_str(Filter$yacas_cmd)
    z <- Ryacas::ysym(z_res)
    stopifnot(z$is_mat)
    Filtert <- as.numeric(Ryacas::yac_str(paste0("Length(Transpose(", z, "))")))
    if (isFALSE(identical(At, Filtert))) {
      stop(
        "`A` and `Filter` do not have compatible dimensions."
      )
    }
    expr <- paste0(z, "*", C, "*", "Transpose(", z, ")")
  } else {
    expr <- C
  }
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
