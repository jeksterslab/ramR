#' Scaled/Standardized RAM Matrices
#'
#' Derives the scaled/standardized RAM matrices.
#'
#' The scaled/standardized \eqn{\mathbf{A}} and \eqn{\mathbf{S}}
#' are  given by
#'
#' \deqn{
#'   \mathbf{A}_{\mathrm{scaled}}
#'   =
#'   \mathbf{D} \mathbf{A} \mathbf{D}^{-1}
#' }
#'
#' \deqn{
#'   \mathbf{S}_{\mathrm{scaled}}
#'   =
#'   \mathbf{D} \mathbf{S} \mathbf{D}
#' }
#'
#' where \eqn{\mathbf{D}} is a diagonal matrix
#' whose diagonal elements are the diagonal elements of \eqn{\mathbf{C}}
#' raised to \eqn{-\frac{1}{2}},
#' that is, the inverse of the standard deviations of the variables.
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a list with the following elements
#'
#'   \describe{
#'     \item{A.scaled}{
#'       `t by t` matrix \eqn{\mathbf{A}_{\mathrm{scaled}}}.
#'       Scaled/standardized asymmetric paths (single-headed arrows),
#'       such as regression coefficients and factor loadings.
#'     }
#'     \item{S.scaled}{
#'       `t by t` numeric matrix \eqn{\mathbf{S}_{\mathrm{scaled}}}.
#'       Scaled/standardized symmetric paths (double-headed arrows),
#'       such as variances and covariances.
#'     }
#'   }
#'
#' @inheritParams CheckRAMMatrices
#' @inheritParams IminusA
#' @export
RAMScaled <- function(A,
                      S,
                      Filter,
                      C = NULL,
                      C.scaled = NULL,
                      check = TRUE,
                      ...) {
  UseMethod("RAMScaled")
}

#' @rdname RAMScaled
#' @inheritParams IminusA
#' @inheritParams RAMScaled
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
#' (RAM <- RAMScaled(A, S, Filter))
#' C(A = RAM$A.scaled, S = RAM$S.scaled)
#' M(A = RAM$A.scaled, S = RAM$S.scaled, Filter = Filter)
#' @export
RAMScaled.default <- function(A,
                              S,
                              Filter,
                              C = NULL,
                              C.scaled = NULL,
                              check = TRUE,
                              ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S,
      Filter = Filter,
      C = C,
      C.scaled = C.scaled
    )
    A <- RAM$A
    S <- RAM$S
    Filter <- RAM$Filter
    C <- RAM$C
    C.scaled <- RAM$C.scaled
  }
  if (is.null(C) || is.null(C.scaled)) {
    Expectations <- Expectations(
      A = A,
      S = S,
      Filter = Filter,
      check = FALSE
    )
    C <- Expectations$C
    C.scaled <- Expectations$C.scaled
  }
  t <- dim(A)[1]
  InvD <- matrix(
    0,
    nrow = t,
    ncol = t
  )
  diag(InvD) <- sqrt(diag(C))
  InvD <- solve(InvD)
  A.scaled <- InvD %*% A %*% solve(InvD)
  S.scaled <- InvD %*% S %*% InvD
  rownames(A.scaled) <- colnames(A.scaled) <- colnames(A)
  rownames(S.scaled) <- colnames(S.scaled) <- colnames(A)
  return(
    list(
      A.scaled = A.scaled,
      S.scaled = S.scaled
    )
  )
}

#' @rdname RAMScaled
#' @inheritParams IminusA
#' @inheritParams RAMScaled
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
#' (RAM <- RAMScaled(Ryacas::ysym(A), S, Filter))
#' RAMScaled(Ryacas::ysym(A), S, Filter, format = "str")
#' RAMScaled(Ryacas::ysym(A), S, Filter, format = "tex")
#' RAMScaled(Ryacas::ysym(A), S, Filter, R = TRUE)
#'
#' C(A = RAM$A.scaled, S = RAM$S.scaled)
#' M(A = RAM$A.scaled, S = RAM$S.scaled, Filter = Filter)
#'
#' # Assigning values to symbols
#'
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#'
#' RAMScaled(Ryacas::ysym(A), S, Filter)
#' RAMScaled(Ryacas::ysym(A), S, Filter, format = "str")
#' RAMScaled(Ryacas::ysym(A), S, Filter, format = "tex")
#' RAMScaled(Ryacas::ysym(A), S, Filter, R = TRUE)
#' eval(RAMScaled(Ryacas::ysym(A), S, Filter, R = TRUE))
#'
#' C(A = RAM$A.scaled, S = RAM$S.scaled)
#' M(A = RAM$A.scaled, S = RAM$S.scaled, Filter = Filter)
#' @export
RAMScaled.yac_symbol <- function(A,
                                 S,
                                 Filter,
                                 C = NULL,
                                 C.scaled = NULL,
                                 check = TRUE,
                                 exe = TRUE,
                                 R = FALSE,
                                 format = "ysym",
                                 simplify = FALSE,
                                 ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S,
      Filter = Filter,
      C = C,
      C.scaled = C.scaled
    )
    A <- RAM$A
    S <- RAM$S
    Filter <- RAM$Filter
    C <- RAM$C
    C.scaled <- RAM$C.scaled
  } else {
    S <- yacR::as.ysym(S)
    if (!is.null(Filter)) {
      Filter <- yacR::as.ysym(Filter)
    }
    if (!is.null(C)) {
      C <- yacR::as.ysym(C)
    }
    if (!is.null(C.scaled)) {
      C.scaled <- yacR::as.ysym(C.scaled)
    }
  }
  if (is.null(C) || is.null(C.scaled)) {
    Expectations <- Expectations(
      A = A,
      S = S,
      Filter = Filter,
      check = FALSE
    )
    C <- Expectations$C
    C.scaled <- Expectations$C.scaled
  }
  InvD <- paste0(
    "Inverse(",
    "DiagonalMatrix(",
    "Sqrt(",
    "Diagonal(",
    C,
    ")",
    ")",
    ")",
    ")"
  )
  A.scaled <- paste0(
    InvD,
    "*",
    A,
    "*",
    "Inverse(",
    InvD,
    ")"
  )
  S.scaled <- paste0(
    InvD,
    "*",
    S,
    "*",
    InvD
  )
  A.scaled <- yacR::as.ysym(A.scaled)
  S.scaled <- yacR::as.ysym(S.scaled)
  if (exe) {
    A.scaled <- yacR::Exe(
      A.scaled,
      R = R,
      format = format,
      simplify = simplify
    )
    S.scaled <- yacR::Exe(
      S.scaled,
      R = R,
      format = format,
      simplify = simplify
    )
  }
  return(
    list(
      A.scaled = A.scaled,
      S.scaled = S.scaled
    )
  )
}
