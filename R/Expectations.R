#' Expectations from the Reticular Action Model (RAM) Matrices
#'
#' Derives the mean and covariance expectations
#' from the Reticular Action Model (RAM) matrices.
#'
#' @return Returns list with the following elements
#'
#'   \describe{
#'     \item{A}{
#'       A `t by t` matrix \eqn{\mathbf{A}}.
#'       Asymmetric paths (single-headed arrows),
#'       such as regression coefficients and factor loadings.
#'     }
#'     \item{S}{
#'       S `t by t` numeric matrix \eqn{\mathbf{S}}.
#'       Symmetric paths (double-headed arrows),
#'       such as variances and covariances.
#'     }
#'     \item{u}{`t by 1` matrix of mean structure parameters.}
#'     \item{Filter}{
#'       Filter `p by t` numeric matrix
#'       \eqn{\mathbf{F}}.
#'       Filter matrix used to select observed variables.
#'     }
#'     \item{v}{`t by 1` matrix of expected values.}
#'     \item{g}{`p by 1` matrix of expected values of observed variables.}
#'     \item{C}{`t by t` matrix of expected covariances.}
#'     \item{M}{`p by p` matrix of expected covariances of observed variables.}
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inherit u details
#' @inherit v details
#' @inherit g details
#' @inherit C details
#' @inherit M details
#'
#' @inheritParams M
#' @inheritParams g
#'
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
#' u <- c(0.00, 0.50, 0.00)
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(Filter) <- c("y", "x", "e")
#' Expectations(A, S, u, Filter)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigmax2", "sigmae2")
#' u <- c("alpha", "mux", 0)
#' Expectations(Ryacas::ysym(A), S, u, Filter)
#' Expectations(Ryacas::ysym(A), S, u, Filter, tex = TRUE)
#' Expectations(Ryacas::ysym(A), S, u, Filter, ysym = FALSE)
#' Expectations(Ryacas::ysym(A), S, u, Filter, str = FALSE)
#'
#' alpha <- 0
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' mux <- 0.50
#' Expectations(Ryacas::ysym(A), S, u, Filter)
#' Expectations(Ryacas::ysym(A), S, u, Filter, tex = TRUE)
#' Expectations(Ryacas::ysym(A), S, u, Filter, ysym = FALSE)
#' Expectations(Ryacas::ysym(A), S, u, Filter, str = FALSE)
#' Expectations <- Expectations(
#'   Ryacas::ysym(A), S, u, Filter,
#'   str = FALSE
#' )
#' eval(Expectations$A)
#' eval(Expectations$S)
#' eval(Expectations$u)
#' eval(Expectations$Filter)
#' eval(Expectations$v)
#' eval(Expectations$g)
#' eval(Expectations$C)
#' eval(Expectations$M)
#' @export
Expectations <- function(A,
                         S,
                         u = NULL,
                         Filter = NULL,
                         ...) {
  UseMethod("Expectations")
}

#' @rdname Expectations
#' @inheritParams IminusA
#' @inheritParams Expectations
#' @export
Expectations.default <- function(A,
                                 S,
                                 u = NULL,
                                 Filter = NULL,
                                 str = TRUE,
                                 ysym = TRUE,
                                 simplify = FALSE,
                                 tex = FALSE,
                                 ...) {
  C <- C.default(
    A = A,
    S = S
  )
  if (is.null(Filter)) {
    M <- C
  } else {
    M <- Filter %*% C %*% t(Filter)
  }
  if (isFALSE(is.null(u))) {
    u <- matrix(
      u,
      ncol = 1
    )
    colnames(u) <- "u"
    rownames(u) <- rownames(A)
    v <- v.default(
      A = A,
      u = u
    )
    if (is.null(Filter)) {
      g <- v
    } else {
      g <- Filter %*% v
    }
    colnames(g) <- "g"
  } else {
    v <- NULL
    g <- NULL
  }
  if (is.null(Filter)) {
    Filter <- matrixR::IdentityFrom(
      A
    )
  }
  return(
    list(
      A = A,
      S = S,
      u = u,
      Filter = Filter,
      v = v,
      g = g,
      C = C,
      M = M
    )
  )
}

#' @rdname Expectations
#' @inheritParams IminusA
#' @inheritParams Expectations
#' @export
Expectations.yac_symbol <- function(A,
                                    S,
                                    u = NULL,
                                    Filter = NULL,
                                    str = TRUE,
                                    ysym = TRUE,
                                    simplify = FALSE,
                                    tex = FALSE,
                                    exe = TRUE,
                                    ...) {
  Cysym <- C(
    A = A,
    S = S,
    exe = FALSE
  )
  Cout <- YacExe(
    expr = Cysym,
    str = str,
    ysym = ysym,
    tex = tex,
    simplify = simplify
  )
  if (is.null(Filter)) {
    Mysym <- Cysym
    Mout <- Cout
  } else {
    Mysym <- M(
      A = A,
      S = S,
      Filter = Filter,
      exe = FALSE
    )
    Mout <- YacExe(
      expr = Mysym,
      str = str,
      ysym = ysym,
      tex = tex,
      simplify = simplify
    )
  }
  if (is.null(u)) {
    vout <- NULL
    gout <- NULL
  } else {
    vysym <- v(
      A = A,
      u = u,
      exe = FALSE
    )
    vout <- YacExe(
      expr = vysym,
      str = str,
      ysym = ysym,
      tex = tex,
      simplify = simplify
    )
    if (is.null(Filter)) {
      gysym <- vysym
      gout <- vout
    } else {
      gysym <- g(
        A = A,
        u = u,
        Filter = Filter,
        exe = FALSE
      )
      gout <- YacExe(
        expr = gysym,
        str = str,
        ysym = ysym,
        tex = tex,
        simplify = simplify
      )
    }
  }
  # make input the same format as output
  Aout <- YacExe(
    expr = A,
    str = str,
    ysym = ysym,
    tex = tex,
    simplify = simplify
  )
  Sout <- YacExe(
    expr = R2Yac(S),
    str = str,
    ysym = ysym,
    tex = tex,
    simplify = simplify
  )
  if (is.null(Filter)) {
    # Filter as identity matrix
    Filter <- diag(
      as.numeric(
        Ryacas::yac_str(
          paste0(
            "Length(",
            A,
            ")"
          )
        )
      )
    )
  }
  Filterout <- YacExe(
    expr = R2Yac(Filter),
    str = str,
    ysym = ysym,
    tex = tex,
    simplify = simplify
  )
  if (is.null(u)) {
    uout <- u
  } else {
    uout <- YacExe(
      expr = RVector2Yac(u),
      str = str,
      ysym = ysym,
      tex = tex,
      simplify = simplify
    )
  }
  return(
    list(
      A = Aout,
      S = Sout,
      u = uout,
      Filter = Filterout,
      v = vout,
      g = gout,
      C = Cout,
      M = Mout
    )
  )
}
