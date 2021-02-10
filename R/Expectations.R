#' Expectations from the Reticular Action Model (RAM) Matrices
#'
#' Derives the means and covariance expectations
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
  if (isFALSE(is.null(Filter))) {
    FilterNull <- TRUE
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
  }
  I <- paste0("Identity(Length(", y, "))")
  E <- paste0("Inverse(", I, "-", y, ")")
  C <- paste0(E, "*", x, "*", "Transpose(", E, ")")
  if (FilterNull) {
    M <- paste0(z, "*", C, "*", "Transpose(", z, ")")
  } else {
    M <- C
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
    a <- Ryacas::ysym(u)
    v <- paste0(E, "*", a)
    if (FilterNull) {
      g <- paste0(Filter, "*", v)
    } else {
      g <- v
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
      A <- y
      S <- x
      Filter <- z
      u <- a
    } else {
      A <- Ryacas::yac_str(y)
      S <- Ryacas::yac_str(x)
      Filter <- Ryacas::yac_str(z)
      u <- Ryacas::yac_str(a)
    }
    if (tex) {
      A <- Ryacas::tex(y)
      S <- Ryacas::tex(x)
      Filter <- Ryacas::tex(z)
      u <- Ryacas::tex(a)
    }
  } else {
    A <- Ryacas::yac_expr(y)
    S <- Ryacas::yac_expr(x)
    Filter <- Ryacas::yac_expr(z)
    u <- Ryacas::yac_expr(a)
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
