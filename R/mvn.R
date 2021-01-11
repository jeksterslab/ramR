#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Generate Data from a Multivariate Normal Distribution
#'   Using the Reticular Action Model (RAM) Notation
#'
#' @description Generates data from a multivariate normal distribution
#'   using the Reticular Action Model (RAM) notation
#'   with the [MASS::mvrnorm()] function.
#'
#' @details The Multivariate Normal Distribution is given by
#'   \deqn{
#'     \mathbf{X} \sim \mathcal{N}_{k} \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right)
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
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \mathbf{M}
#'   }
#'   \deqn{
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \boldsymbol{\Omega}
#'     \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}}
#'   }
#'
#' @inheritParams Sigmatheta
#' @inheritParams mutheta
#' @param n Integer.
#'   Sample size.
#' @param ... Additional arguments to pass to [MASS::mvrnorm()].
#' @export
mvn <- function(n,
                A,
                Omega,
                M,
                filter = NULL,
                ...) {
  mu <- mutheta(
    M = M,
    A = A,
    filter = filter
  )
  Sigma <- Sigmatheta(
    A = A,
    Omega = Omega,
    filter = filter
  )
  return(
    MASS::mvrnorm(
      n = n,
      mu = mu,
      Sigma = Sigma,
      ...
    )
  )
}
