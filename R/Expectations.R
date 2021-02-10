#' Expectations from the Reticular Action Model (RAM) Matrices
#'
#' Derives the mean and covariance expectations
#' from the Reticular Action Model (RAM) matrices.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inherit u details
#' @inherit v details
#' @inherit g details
#' @inherit C details
#' @inherit M details
#' @inheritParams M
#' @inheritParams g
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
#' Expectations <- Expectations(Ryacas::ysym(A), S, u, Filter, str = FALSE)
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
    Filter <- matrixR::IdentityFrom(A)
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
                                    ...) {
  stopifnot(methods::is(A, "yac_symbol"))
  Aysym <- Ryacas::ysym(Ryacas::yac_str(A$yacas_cmd))
  stopifnot(Aysym$is_mat)
  stopifnot(matrixR::IsSquareMatrix(Aysym))
  # apply IsNilpotent in the future
  if (methods::is(S, "yac_symbol")) {
    Sysym <- S
  } else {
    Sysym <- Ryacas::ysym(S)
  }
  Sysym <- Ryacas::ysym(Ryacas::yac_str(Sysym$yacas_cmd))
  stopifnot(Sysym$is_mat)
  stopifnot(matrixR::IsSymmetric(Sysym))
  ADimensions <- as.numeric(Ryacas::yac_str(paste0("Length(", Aysym, ")")))
  SDimensions <- as.numeric(Ryacas::yac_str(paste0("Length(", Sysym, ")")))
  stopifnot(identical(ADimensions, SDimensions))
  I <- paste0("Identity(Length(", Aysym, "))")
  E <- paste0("Inverse(", I, "-", Aysym, ")")
  C <- paste0(E, "*", Sysym, "*", "Transpose(", E, ")")
  if (is.null(Filter)) {
    M <- C
    Filterysym <- Ryacas::ysym(diag(ADimensions))
  } else {
    if (methods::is(Filter, "yac_symbol")) {
      Filterysym <- Filter
    } else {
      Filterysym <- Ryacas::ysym(Filter)
    }
    Filterysym <- Ryacas::ysym(Ryacas::yac_str(Filterysym$yacas_cmd))
    stopifnot(Filterysym$is_mat)
    FilterDimensions <- as.numeric(Ryacas::yac_str(paste0("Length(Transpose(", Filterysym, "))")))
    stopifnot(identical(ADimensions, FilterDimensions))
    M <- paste0(Filterysym, "*", C, "*", "Transpose(", Filterysym, ")")
  }
  C <- .exe(
    expr = C,
    str = str,
    ysym = ysym,
    simplify = simplify,
    tex = tex
  )
  M <- .exe(
    expr = M,
    str = str,
    ysym = ysym,
    simplify = simplify,
    tex = tex
  )
  if (is.null(u)) {
    v <- NULL
    g <- NULL
  } else {
    u <- matrix(
      u,
      ncol = 1
    )
    uysym <- Ryacas::ysym(u)
    v <- paste0(E, "*", uysym)
    if (is.null(Filter)) {
      g <- v
    } else {
      g <- paste0(Filterysym, "*", v)
    }
    v <- .exe(
      expr = v,
      str = str,
      ysym = ysym,
      simplify = simplify,
      tex = tex
    )
    g <- .exe(
      expr = g,
      str = str,
      ysym = ysym,
      simplify = simplify,
      tex = tex
    )
  }
  if (str) {
    if (ysym) {
      A <- Aysym
      S <- Sysym
      Filter <- Filterysym
      if (isFALSE(is.null(u))) {
        u <- uysym
      }
    } else {
      A <- Ryacas::yac_str(Aysym)
      S <- Ryacas::yac_str(Sysym)
      Filter <- Ryacas::yac_str(Filterysym)
      if (isFALSE(is.null(u))) {
        u <- Ryacas::yac_str(uysym)
      }
    }
    if (tex) {
      A <- Ryacas::tex(Aysym)
      S <- Ryacas::tex(Sysym)
      Filter <- Ryacas::tex(Filterysym)
      if (isFALSE(is.null(u))) {
        u <- Ryacas::tex(uysym)
      }
    }
  } else {
    A <- Ryacas::yac_expr(Aysym)
    S <- Ryacas::yac_expr(Sysym)
    Filter <- Ryacas::yac_expr(Filterysym)
    if (isFALSE(is.null(u))) {
      u <- Ryacas::yac_expr(uysym)
    }
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
