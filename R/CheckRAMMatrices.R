#' Check Matrices
#'
#' Checks if the input matrices are in the correct form.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family preprocess functions
#' @keywords preprocess
#'
#' @param A `t by t` matrix \eqn{\mathbf{A}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param S `t by t` numeric matrix \eqn{\mathbf{S}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @param u vector of length `t` or `t by 1` matrix.
#'   Mean structure parameters.
#' @param Filter `p by t` numeric matrix
#'   \eqn{\mathbf{F}}.
#'   Filter matrix used to select observed variables.
#' @param C `t by t` numeric matrix \eqn{\mathbf{C}}.
#'   Model-implied variance-covariance matrix.
#' @param C.scaled `t by t` numeric matrix
#'   \eqn{\mathbf{C}_{\mathrm{scaled}}}.
#'   Scaled/standardized model-implied variance-covariance matrix.
#' @param v vector of length `t` or `t by 1` matrix.
#'   Expected values.
#' @param ... ...
#' @export
CheckRAMMatrices <- function(A,
                             S = NULL,
                             u = NULL,
                             Filter = NULL,
                             C = NULL,
                             C.scaled = NULL,
                             v = NULL,
                             ...) {
  UseMethod("CheckRAMMatrices")
}

#' @rdname CheckRAMMatrices
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' diag(S) <- c(0, 0.25, 1)
#' u <- c(0.00, 0.50, 0.00)
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' C <- matrix(
#'   data = c(
#'     1.25, 0.25, 1.00,
#'     0.25, 0.25, 0.00,
#'     1.00, 0.00, 1.00
#'   ),
#'   nrow = dim(A)[1]
#' )
#' v <- c(0.50, 0.50, 0.00)
#' CheckRAMMatrices(
#'   A = A,
#'   S = S,
#'   u = u,
#'   Filter = Filter,
#'   C = C,
#'   v = v
#' )
#' @export
CheckRAMMatrices.default <- function(A,
                                     S = NULL,
                                     u = NULL,
                                     Filter = NULL,
                                     C = NULL,
                                     C.scaled = NULL,
                                     v = NULL,
                                     ...) {
  stopifnot(is.numeric(A))
  stopifnot(matrixR::IsNilpotent(A))
  A_dimensions <- dim(A)
  names <- colnames(A)
  if (!is.null(S)) {
    stopifnot(is.numeric(S))
    stopifnot(matrixR::IsSymmetric(round(S, digits = 4)))
    S_dimensions <- dim(S)
    if (!is.null(names)) {
      colnames(S) <- rownames(S) <- names
    }
  }
  if (!is.null(u)) {
    stopifnot(is.numeric(u))
    u <- as.matrix(u)
    u_dimensions <- dim(u)
    colnames(u) <- "u"
    if (!is.null(names)) {
      rownames(u) <- names
    }
  }
  if (!is.null(Filter)) {
    stopifnot(is.numeric(Filter))
    Filter_dimensions <- dim(Filter)
    if (!is.null(names)) {
      colnames(Filter) <- names
      rows <- colSums(Filter)
      rownames(Filter) <- names(rows[rows == 1])
    }
  }
  if (!is.null(C)) {
    stopifnot(is.numeric(C))
    stopifnot(matrixR::IsSymmetric(round(C, digits = 4)))
    C_dimensions <- dim(C)
    if (!is.null(names)) {
      colnames(C) <- rownames(C) <- names
    }
  }
  if (!is.null(C.scaled)) {
    stopifnot(is.numeric(C.scaled))
    stopifnot(matrixR::IsSymmetric(round(C.scaled, digits = 4)))
    C.scaled_dimensions <- dim(C.scaled)
    if (!is.null(names)) {
      colnames(C.scaled) <- rownames(C.scaled) <- names
    }
  }
  if (!is.null(v)) {
    stopifnot(is.numeric(v))
    v <- as.matrix(v)
    v_dimensions <- dim(v)
    colnames(v) <- "v"
    if (!is.null(names)) {
      rownames(v) <- names
    }
  }
  if (!is.null(S)) {
    stopifnot(
      identical(
        S_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(u)) {
    stopifnot(
      identical(
        u_dimensions[1],
        A_dimensions[1]
      )
    )
  }
  if (!is.null(Filter)) {
    stopifnot(
      identical(
        Filter_dimensions[2],
        A_dimensions[2]
      )
    )
  }
  if (!is.null(C)) {
    stopifnot(
      identical(
        C_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(C.scaled)) {
    stopifnot(
      identical(
        C.scaled_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(v)) {
    stopifnot(
      identical(
        v_dimensions[1],
        A_dimensions[1]
      )
    )
  }
  return(
    list(
      A = A,
      S = S,
      u = u,
      Filter = Filter,
      C = C,
      C.scaled = C.scaled,
      v = v
    )
  )
}

#' @rdname CheckRAMMatrices
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigmax2", "sigmae2")
#' u <- c("alpha", "mux", 0)
#' C <- matrix(
#'   data = c(
#'     "sigmax2*beta^2+sigmae2", "sigmax2*beta", "sigmae2",
#'     "sigmax2*beta", "sigmax2", 0,
#'     "sigmae2", 0, "sigmae2"
#'   ),
#'   nrow = dim(A)[1]
#' )
#' v <- c("alpha+beta*mux", "mux", 0)
#' CheckRAMMatrices(
#'   A = Ryacas::ysym(A),
#'   S = S,
#'   u = u,
#'   Filter = Filter,
#'   C = C,
#'   v = v
#' )
#' @export
CheckRAMMatrices.yac_symbol <- function(A,
                                        S = NULL,
                                        u = NULL,
                                        Filter = NULL,
                                        C = NULL,
                                        C.scaled = NULL,
                                        v = NULL,
                                        ...) {
  stopifnot(matrixR::IsSquareMatrix(A))
  A_dimensions <- Length(A)
  if (!is.null(S)) {
    # IsSymmetric is not safe for symbolic results.
    # See https://github.com/grzegorzmazur/yacas/issues/327
    stopifnot(matrixR::IsSquareMatrix(S))
    S <- yacR::as.ysym.mat(S)
    S_dimensions <- Length(S)
  }
  if (!is.null(u)) {
    u <- yacR::as.ysym.mat(u)
    u_dimensions <- Length(u)
  }
  if (!is.null(Filter)) {
    Filter <- yacR::as.ysym.mat(Filter)
    Filter_dimensions <- Length(Filter, transpose = TRUE)
  }
  if (!is.null(C)) {
    # IsSymmetric is not safe for symbolic results.
    # See https://github.com/grzegorzmazur/yacas/issues/327
    stopifnot(matrixR::IsSquareMatrix(C))
    C <- yacR::as.ysym.mat(C)
    C_dimensions <- Length(C)
  }
  if (!is.null(C.scaled)) {
    # IsSymmetric is not safe for symbolic results.
    # See https://github.com/grzegorzmazur/yacas/issues/327
    stopifnot(matrixR::IsSquareMatrix(C.scaled))
    C.scaled <- yacR::as.ysym.mat(C.scaled)
    C.scaled_dimensions <- Length(C.scaled)
  }
  if (!is.null(v)) {
    v <- yacR::as.ysym.mat(v)
    v_dimensions <- Length(v)
  }
  if (!is.null(S)) {
    stopifnot(
      identical(
        S_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(u)) {
    stopifnot(
      identical(
        u_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(Filter)) {
    stopifnot(
      identical(
        Filter_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(C)) {
    stopifnot(
      identical(
        C_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(C.scaled)) {
    stopifnot(
      identical(
        C.scaled_dimensions,
        A_dimensions
      )
    )
  }
  if (!is.null(v)) {
    stopifnot(
      identical(
        v_dimensions,
        A_dimensions
      )
    )
  }
  return(
    list(
      A = A,
      S = S,
      u = u,
      Filter = Filter,
      C = C,
      C.scaled = C.scaled,
      v = v
    )
  )
}
