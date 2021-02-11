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
#' @author Ivan Jacob Agaloos Pesigan
#' @family RAM matrices functions
#' @keywords ram
#' @inherit ramR references
#' @inheritParams v
#' @inheritParams M
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' # Numeric -----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' u <- c(0.00, 0.50, 0.00)
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(A) <- rownames(A) <- c("y", "x", "e")
#' g(A, u, Filter)
#'
#' # Symbolic ----------------------------------------------------------
#' A <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, "beta", 1)
#' u <- c("alpha", "mux", 0)
#' g(Ryacas::ysym(A), u, Filter)
#' g(Ryacas::ysym(A), u, Filter, tex = TRUE)
#' g(Ryacas::ysym(A), u, Filter, ysym = FALSE)
#' g(Ryacas::ysym(A), u, Filter, str = FALSE)
#'
#' alpha <- 0
#' beta <- 1
#' mux <- 0.50
#' g(Ryacas::ysym(A), u, Filter)
#' g(Ryacas::ysym(A), u, Filter, tex = TRUE)
#' g(Ryacas::ysym(A), u, Filter, ysym = FALSE)
#' g(Ryacas::ysym(A), u, Filter, str = FALSE)
#' eval(g(Ryacas::ysym(A), u, Filter, str = FALSE))
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
#' @export
g.default <- function(A,
                      u,
                      Filter = NULL,
                      str = TRUE,
                      ysym = TRUE,
                      simplify = FALSE,
                      tex = FALSE,
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
#' @export
g.yac_symbol <- function(A,
                         u,
                         Filter = NULL,
                         str = TRUE,
                         ysym = TRUE,
                         simplify = FALSE,
                         tex = FALSE,
                         ...) {
  stopifnot(
    methods::is(
      A,
      "yac_symbol"
    )
  )
  Aysym <- Ryacas::ysym(
    Ryacas::yac_str(
      A$yacas_cmd
    )
  )
  stopifnot(
    Aysym$is_mat
  )
  stopifnot(
    matrixR::IsSquareMatrix(
      Aysym
    )
  )
  I <- paste0(
    "Identity(Length(",
    Aysym,
    "))"
  )
  E <- paste0(
    "Inverse(",
    I,
    "-",
    Aysym,
    ")"
  )
  v <- paste0(
    E,
    "*",
    Ryacas::ysym(
      matrix(
        u,
        ncol = 1
      )
    )
  )
  if (is.null(Filter)) {
    expr <- v
  } else {
    if (methods::is(Filter, "yac_symbol")) {
      Filterysym <- Filter
    } else {
      Filterysym <- Ryacas::ysym(
        Filter
      )
    }
    Filterysym <- Ryacas::ysym(
      Ryacas::yac_str(
        Filterysym$yacas_cmd
      )
    )
    stopifnot(
      Filterysym$is_mat
    )
    ADimensions <- as.numeric(
      Ryacas::yac_str(
        paste0(
          "Length(",
          Aysym,
          ")"
        )
      )
    )
    FilterDimensions <- as.numeric(
      Ryacas::yac_str(
        paste0(
          "Length(Transpose(",
          Filterysym,
          "))"
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
      v
    )
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
