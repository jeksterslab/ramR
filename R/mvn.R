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
#' @inheritParams M_num
#' @inheritParams g_num
#' @param n Integer.
#'   Sample size.
#' @param ... Additional arguments to pass to [MASS::mvrnorm()].
#' @examples
#' # This is an example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' n <- 1000
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, 2] <- 1
#' A[1, 3] <- 1
#' diag(S) <- c(0, 0.25, 0.25)
#' u <- c(-0.50, 0.50, 0.00)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' colnames(filter) <- c("y", "x", "e")
#' rownames(filter) <- c("y", "x")
#' mvn(n, A, S, u, filter)
#' @export
mvn <- function(n,
                A,
                S,
                u,
                filter,
                ...) {
  v <- v_num(
    A,
    u
  )
  g <- filter %*% v
  colnames(g) <- "g"
  C <- C_num(
    A,
    S
  )
  M <- filter %*% C %*% t(filter)
  out <- MASS::mvrnorm(
    n = n,
    mu = g,
    Sigma = M,
    ...
  )
  attributes(out)$n <- n
  attributes(out)$A <- A
  attributes(out)$S <- S
  attributes(out)$u <- u
  attributes(out)$filter <- filter
  attributes(out)$v <- v
  attributes(out)$g <- g
  attributes(out)$C <- C
  attributes(out)$M <- M
  return(
    out
  )
}
