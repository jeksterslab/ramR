#' Scaled/Standardized RAM Matrices
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
#' @param C.scaled `t by t` numeric matrix
#'   \eqn{\mathbf{C}_{\mathrm{scaled}}}.
#'   Scaled/standardized model-implied variance-covariance matrix.
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
      C = C
    )
    A <- RAM$A
    S <- RAM$S
    Filter <- RAM$Filter
    C <- RAM$C
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
  # t <- dim(A)[1]
  # S.scaled <- A.scaled <- matrix(
  #  0,
  #  nrow = t,
  #  ncol = t
  # )
  # colnames(A.scaled) <- rownames(A.scaled) <- colnames(A)
  # colnames(S.scaled) <- rownames(S.scaled) <- colnames(A)
  # for (x in 1:t) {
  #  for (y in 1:t) {
  #    if (A[y, x] != 0) {
  #      A.scaled[y, x] <- A[y, x] * sqrt(C[x, x]) / sqrt(C[y, y])
  #    }
  #  }
  # }
  ## diagonal elements
  # for (y in 1:t) {
  #  if (S[y, y] != 0) {
  #    S.scaled[y, y] <- 1 - sum(A.scaled[y, ]^2)
  #  }
  # }
  ## off-diagonal elements
  # for (x in 1:t) {
  #  for (y in 1:t) {
  #    if (x != y) {
  #      if (S[y, x] != 0) {
  #        S.scaled[y, x] <- C.scaled[y, x]
  #      }
  #    }
  #  }
  # }
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
      C = C
    )
    A <- RAM$A
    S <- RAM$S
    Filter <- RAM$Filter
    C <- RAM$C
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
  # t <- Length(A)
  # S.scaled <- A.scaled <- matrix(
  #  0,
  #  nrow = t,
  #  ncol = t
  # )
  # for (x in 1:t) {
  #  for (y in 1:t) {
  #    a <- suppressWarnings(as.numeric(A[y, x]$yacas_cmd))
  #    if (is.na(a) || a != 0) {
  #      A.scaled[y, x] <- paste0(
  #        "(",
  #        "(",
  #        A[y, x]$yacas_cmd,
  #        ")",
  #        "*",
  #        "Sqrt(",
  #        C[x, x]$yacas_cmd,
  #        ")",
  #        "/",
  #        "Sqrt(",
  #        C[y, y]$yacas_cmd,
  #        ")",
  #        ")"
  #      )
  #    }
  #  }
  # }
  ## diagonal elements
  # for (y in 1:t) {
  #  s <- suppressWarnings(as.numeric(S[y, y]$yacas_cmd))
  #  if (is.na(s) || s != 0) {
  #    S.scaled[y, y] <- paste0(
  #      "1",
  #      "-",
  #      "(",
  #      paste0(
  #        paste0(
  #          A.scaled[y, ],
  #          "^2"
  #        ),
  #        collapse = "+"
  #      ),
  #      ")"
  #    )
  #  }
  # }
  ## off-diagonal elements
  # for (x in 1:t) {
  #  for (y in 1:t) {
  #    if (x != y) {
  #      s <- suppressWarnings(as.numeric(S[y, x]$yacas_cmd))
  #      if (is.na(s) || s != 0) {
  #        S.scaled[y, x] <- C.scaled[y, x]$yacas_cmd
  #      }
  #    }
  #  }
  # }
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
