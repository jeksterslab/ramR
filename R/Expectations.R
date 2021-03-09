#' Expectations from the Reticular Action Model (RAM) Matrices
#'
#' Derives the mean and covariance expectations
#' from the Reticular Action Model (RAM) matrices.
#'
#' The vector of expected values \eqn{\mathbf{v}}
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#'
#'   \deqn{
#'     \mathbf{v}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}}
#'     \mathbf{u} \\
#'     =
#'     \mathbf{E}
#'     \mathbf{u}
#'   }
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
#' The matrix of covariance expectations \eqn{\mathbf{C}}
#' as a function of Reticular Action Model (RAM) matrices is given by
#'
#'   \deqn{
#'     \mathbf{C}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\mathsf{T}} \\
#'     =
#'     \mathbf{E}
#'     \mathbf{S}
#'     \mathbf{E}^{\mathsf{T}}
#'   }
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
#' The matrix of scaled/standardized covariance expectations
#' \eqn{\mathbf{C}_{\mathrm{scaled}}}
#' (also known as correlations)
#' is given by
#'
#'   \deqn{
#'     \mathbf{C}_{\mathrm{scaled}}
#'     =
#'     \mathbf{D}^{-1}
#'     \mathbf{C}
#'     \mathbf{D}^{-1}
#'   }
#'
#' The matrix of scaled/standardized covariance expectations
#' for given variables
#' \eqn{\mathbf{M}_{\mathrm{scaled}}}
#' (also known as correlations)
#' is given by
#'
#'   \deqn{
#'     \mathbf{M}_{\mathrm{scaled}}
#'     =
#'     \mathbf{F}
#'     \mathbf{C}_{\mathrm{scaled}}
#'     \mathbf{F}^{\mathsf{T}}
#'   }
#'
#' where
#'
#'   - \eqn{\mathbf{u}_{t \times 1}} vector of parameters
#'     for the mean structure,
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{F}_{p \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\mathbf{D}_{t \times t}} represents a diagonal matrix
#'     who diagonal elements are the square root of the diagonal elements
#'     of \eqn{\mathbf{C}},
#'   - \eqn{p} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables,
#'     that is \eqn{p + q} .
#'
#' @return Returns a list with the following elements
#'
#'   \describe{
#'     \item{v}{
#'       `t by 1` matrix \eqn{\mathbf{v}}
#'       of expected values.
#'     }
#'     \item{g}{
#'       `p by 1` matrix \eqn{\mathbf{g}}
#'       of expected values of observed variables.
#'     }
#'     \item{C}{
#'       `t by t` matrix \eqn{\mathbf{C}}
#'       of expected covariances.
#'     }
#'     \item{M}{
#'       `p by p` matrix \eqn{\mathbf{M}}
#'       of expected covariances of observed variables.
#'     }
#'     \item{C.scaled}{
#'       `t by t` matrix \eqn{\mathbf{C}}
#'       of scaled/standardized expected covariances
#'       (also known as correlations).
#'     }
#'     \item{M.scaled}{
#'       `p by p` matrix \eqn{\mathbf{M}}
#'       of scaled/standardized expected covariances
#'       (also known as correlations)
#'       of observed variables.
#'     }
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family RAM matrices functions
#' @keywords ram
#'
#' @inherit ramR references
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
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "ysym")
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "str")
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "tex")
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = TRUE)
#'
#' # Assigning values to symbols
#'
#' alpha <- 0
#' beta <- 1
#' sigmax2 <- 0.25
#' sigmae2 <- 1
#' mux <- 0.50
#'
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "ysym")
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "str")
#' Expectations(Ryacas::ysym(A), S, u, Filter, R = FALSE, format = "tex")
#' (Expectations <- Expectations(Ryacas::ysym(A), S, u, Filter, R = TRUE))
#' eval(Expectations$v)
#' eval(Expectations$g)
#' eval(Expectations$C)
#' eval(Expectations$M)
#' eval(Expectations$C.scaled)
#' eval(Expectations$M.scaled)
#' @export
Expectations <- function(A,
                         S,
                         u = NULL,
                         Filter = NULL,
                         check = TRUE,
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
                                 check = TRUE,
                                 ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S,
      u = u,
      Filter = Filter
    )
    A <- RAM$A
    S <- RAM$S
    u <- RAM$u
    Filter <- RAM$Filter
  }
  C <- C(
    A = A,
    S = S,
    check = FALSE
  )
  C_scaled <- stats::cov2cor(C)
  if (is.null(Filter)) {
    M <- C
    M_scaled <- C_scaled
  } else {
    M <- Filter %*% C %*% t(Filter)
    M_scaled <- Filter %*% C_scaled %*% t(Filter)
  }
  if (!(is.null(u))) {
    v <- v(
      A = A,
      u = u,
      check = FALSE
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
  # if (is.null(Filter)) {
  #  Filter <- matrixR::IdentityFrom(
  #    A
  #  )
  # }
  # RAM.scaled <- RAMScaled(
  #  A = A,
  #  S = S,
  #  Filter = Filter,
  #  C = C,
  #  C.scaled = C_scaled,
  #  check = FALSE
  # )
  return(
    list(
      # A = A,
      # S = S,
      # u = u,
      # Filter = Filter,
      v = v,
      g = g,
      C = C,
      M = M,
      C.scaled = C_scaled,
      M.scaled = M_scaled
      # A.scaled = RAM.scaled$A.scaled,
      # S.scaled = RAM.scaled$S.scaled
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
                                    check = TRUE,
                                    R = FALSE,
                                    format = "ysym",
                                    simplify = FALSE,
                                    ...) {
  if (check) {
    RAM <- CheckRAMMatrices(
      A = A,
      S = S,
      u = u,
      Filter = Filter
    )
    A <- RAM$A
    S <- RAM$S
    u <- RAM$u
    Filter <- RAM$Filter
  } else {
    S <- yacR::as.ysym.mat(S)
    if (!is.null(Filter)) {
      Filter <- yacR::as.ysym.mat(Filter)
    }
    if (!is.null(u)) {
      u <- yacR::as.ysym.mat(u)
    }
  }
  C <- C(
    A = A,
    S = S,
    check = FALSE,
    exe = FALSE
  ) # A and S already checked
  C_out <- yacR::Exe(
    C,
    R = R,
    format = format,
    simplify = simplify
  )
  C_scaled <- yacR::Exe(
    paste0(
      "Inverse(",
      "DiagonalMatrix(",
      "Sqrt(",
      "Diagonal(",
      C,
      ")",
      ")",
      ")",
      ")",
      "*",
      C,
      "*",
      "Inverse(",
      "DiagonalMatrix(",
      "Sqrt(",
      "Diagonal(",
      C,
      ")",
      ")",
      ")",
      ")"
    ),
    R = FALSE,
    format = "ysym",
    simplify = FALSE
  )
  C_scaled_out <- yacR::Exe(
    C_scaled,
    R = R,
    format = format,
    simplify = simplify
  )
  if (is.null(Filter)) {
    M_out <- C_out
    M_scaled_out <- C_scaled_out
  } else {
    M_out <- yacR::Exe(
      paste0(
        Filter,
        "*",
        C,
        "*",
        "Transpose(",
        Filter,
        ")"
      ),
      R = R,
      format = format,
      simplify = simplify
    )
    M_scaled_out <- yacR::Exe(
      paste0(
        Filter,
        "*",
        C_scaled,
        "*",
        "Transpose(",
        Filter,
        ")"
      ),
      R = R,
      format = format,
      simplify = simplify
    )
  }
  if (is.null(u)) {
    v_out <- NULL
    g_out <- NULL
  } else {
    v <- v(
      A = A,
      u = u,
      check = FALSE,
      exe = FALSE
    ) # A and u are already checked
    v_out <- yacR::Exe(
      v,
      R = R,
      format = format,
      simplify = simplify
    )
    if (is.null(Filter)) {
      g <- v
      g_out <- v_out
    } else {
      g <- g(
        A = A,
        u = u,
        Filter = Filter,
        check = FALSE,
        exe = FALSE
      ) # A, u, and Filter are already checked
      g_out <- yacR::Exe(
        g,
        R = R,
        format = format,
        simplify = simplify
      )
    }
  }
  #  # make input the same format as output
  #  A.out <- yacR::Exe(
  #    A,
  #    R = R,
  #    format = format,
  #    simplify = simplify
  #  )
  #  S.out <- yacR::Exe(
  #    yacR::as.ysym(S),
  #    R = R,
  #    format = format,
  #    simplify = simplify
  #  )
  #  if (is.null(Filter)) {
  #    # Filter as identity matrix
  #    Filter <- diag(
  #      as.numeric(
  #        Ryacas::yac_str(
  #          paste0(
  #            "Length(",
  #            A,
  #            ")"
  #          )
  #        )
  #      )
  #    )
  #  }
  #  Filter.out <- yacR::Exe(
  #    yacR::as.ysym(Filter),
  #    R = R,
  #    format = format,
  #    simplify = simplify
  #  )
  #  if (is.null(u)) {
  #    u.out <- u
  #  } else {
  #    u.out <- yacR::Exe(
  #      yacR::as.ysym.mat(u),
  #      R = R,
  #      format = format,
  #      simplify = simplify
  #    )
  #  }
  # RAM.scaled <- RAMScaled(
  #  A = A,
  #  S = S,
  #  Filter = Filter,
  #  C = C,
  #  C.scaled = C_scaled,
  #  check = FALSE,
  #  exe = exe,
  #  R = R,
  #  format = format,
  #  simplify = simplify
  # )
  return(
    list(
      # A = A.out,
      # S = S.out,
      # u = u.out,
      # Filter = Filter.out,
      v = v_out,
      g = g_out,
      C = C_out,
      M = M_out,
      C.scaled = C_scaled_out,
      M.scaled = M_scaled_out
      # A.scaled = RAM.scaled$A.scaled,
      # S.scaled = RAM.scaled$S.scaled
    )
  )
}
