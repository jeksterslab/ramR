#' Vector of Expected Values of Observed Variables \eqn{\mathbf{g}}
#'
#' Derives the vector of expected values
#' of observed variables \eqn{\mathbf{g}}
#' using the Reticular Action Model (RAM) notation.
#'
#' The vector of expected values
#' of observed variables \eqn{\mathbf{g}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{g}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{E}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{F}
#'     \mathbf{v}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\mathbf{u}_{t \times 1}} vector of parameters
#'     for the mean structure,
#'   - \eqn{\mathbf{F}_{p \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{p} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables,
#'     that is \eqn{p + q} .
#'
#' @return \eqn{\mathbf{g} = \mathbf{F} \mathbf{v}}
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
#' @inheritParams v
#' @inheritParams M
#' @export
g <- function(A,
              u,
              Filter = NULL,
              ...) {
  UseMethod("g")
}

#' @rdname g
#' @inheritParams IminusA
#' @inheritParams g
#' @examples
#' # Numeric -----------------------------------------------------------
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' u <- c(0.00, 0.50, 0.00)
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' g(A, u, Filter)
#' @export
g.default <- function(A,
                      u,
                      Filter = NULL,
                      ...) {
  v <- v.default(
    A = A,
    u = u
  )
  if (isFALSE(
    is.null(Filter)
  )
  ) {
    out <- Filter %*% v
    colnames(out) <- "g"
    return(
      out
    )
  } else {
    colnames(v) <- "g"
    return(v)
  }
}

#' @rdname g
#' @inheritParams IminusA
#' @inheritParams g
#' @examples
#' # Symbolic ----------------------------------------------------------
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' u <- c("alpha", "mux", 0)
#' g(Ryacas::ysym(A), u, Filter)
#' g(Ryacas::ysym(A), u, Filter, format = "tex")
#' g(Ryacas::ysym(A), u, Filter, format = "ysym")
#' g(Ryacas::ysym(A), u, Filter, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' g(Ryacas::ysym(A), u, Filter)
#' g(Ryacas::ysym(A), u, Filter, format = "tex")
#' g(Ryacas::ysym(A), u, Filter, format = "ysym")
#' g(Ryacas::ysym(A), u, Filter, R = TRUE)
#' eval(g(Ryacas::ysym(A), u, Filter, R = TRUE))
#' @export
g.yac_symbol <- function(A,
                         u,
                         Filter = NULL,
                         exe = TRUE,
                         R = FALSE,
                         format = "ysym",
                         simplify = FALSE,
                         ...) {
  vysym <- v(
    A = A,
    u = u,
    exe = FALSE
  )
  if (is.null(Filter)) {
    expr <- vysym
  } else {
    Filterysym <- RMatrix2Yac(Filter)
    FilterDimensions <- as.numeric(
      Ryacas::yac_str(
        paste0(
          "Length(Transpose(",
          Filterysym,
          "))"
        )
      )
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
    stopifnot(
      identical(
        ADimensions,
        FilterDimensions
      )
    )
    expr <- paste0(
      Filterysym,
      "*",
      vysym
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
