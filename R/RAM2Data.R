#' Generate Data from a Multivariate Normal Distribution
#'   Using the Reticular Action Model (RAM) Notation
#'
#' Generates data from a multivariate normal distribution
#'   using the Reticular Action Model (RAM) notation
#'   with the [MASS::mvrnorm()] function.
#'
#' The Multivariate Normal Distribution is given by
#'   \deqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left(
#'       \boldsymbol{\mu},
#'       \boldsymbol{\Sigma}
#'     \right)
#'   }
#'
#'   with location parameter
#'
#'   \deqn{
#'     \boldsymbol{\mu} \in \mathbf{R}^{k}
#'   }
#'
#'   and a positive definite covariance matrix
#'
#'   \deqn{
#'     \boldsymbol{\Sigma} \in \mathbf{R}^{k \times k} .
#'   }
#'
#'   The probability density function is given by
#'
#'   \deqn{
#'     f_{\mathbf{X}} \left( x_1, \cdots, x_k \right)
#'     =
#'     \frac{
#'       \exp
#'       \left[
#'         - \frac{1}{2}
#'         \left(
#'           \mathbf{x} - \boldsymbol{\mu}
#'         \right)^{\mathsf{T}}
#'         \boldsymbol{\Sigma}^{-1}
#'         \left(
#'           \mathbf{x} - \boldsymbol{\mu}
#'         \right)
#'       \right]
#'       }{
#'         \sqrt{
#'           \left( 2 \pi \right)^k
#'           | \boldsymbol{\Sigma} |
#'         }
#'       }
#'     }
#'
#'   In this function,
#'   the model-implied mean vector and variance-covariance matrix
#'   are used as parameters to generate the data.
#'
#'   \deqn{
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{g}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{\mathsf{T}}
#'     \mathbf{u}
#'     =
#'     \mathbf{F}
#'     \mathbf{E}
#'     \mathbf{u}
#'   }
#'   \deqn{
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     =
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
#' @author Ivan Jacob Agaloos Pesigan
#' @family data generation functions
#' @keywords data
#' @inherit ramR references
#' @inheritParams Expectations
#' @param n Integer.
#'   Sample size.
#' @param ... Additional arguments to pass to [MASS::mvrnorm()].
#' @examples
#' # This is an example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' n <- 50
#' A <- S <- matrixR::ZeroMatrix(3)
#' A[1, ] <- c(0, 1, 1)
#' diag(S) <- c(0, 0.25, 1)
#' u <- c(0.00, 0.50, 0.00)
#' Filter <- diag(2)
#' Filter <- cbind(Filter, 0)
#' colnames(Filter) <- c("y", "x", "e")
#' rownames(Filter) <- c("y", "x")
#' RAM2Data(n, A, S, u, Filter)
#' @export
RAM2Data <- function(n,
                     A,
                     S,
                     u,
                     Filter,
                     ...) {
  Expectations <- Expectations.default(
    A = A,
    S = S,
    u = u,
    Filter = Filter
  )
  g <- as.vector(
    Expectations$g
  )
  names(g) <- rownames(Filter)
  return(
    MASS::mvrnorm(
      n = n,
      mu = g,
      Sigma = Expectations$M,
      ...
    )
  )
}
