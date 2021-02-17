#' Equations to Multivariate Normal Sample Data
#'
#' Generates data from a multivariate normal distribution
#' from model equations using the [MASS::mvrnorm()] function.
#'
#' The input is a character string
#' that specifies the associations between the variables.
#'
#' The multivariate normal distribution is given by
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
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs <space> op <space> rhs <space> par.label <\n> or <;>`
#'
#'   \describe{
#'     \item{lhs}{
#'       is the variable on the **left-hand side**,
#'     }
#'     \item{rhs}{
#'       is the variable on the **right-hand side**,
#'     }
#'     \item{op}{
#'       is the **operation** between `lhs` and `rhs`,
#'     }
#'     \item{par.label}{
#'       is the column of **parameter label**,
#'     }
#'     \item{\\n or ;}{
#'       are **line breaks**.
#'       **Each line should end with a line break.**
#'     }
#'   }
#'
#' @section Operations:
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{
#'       `left-hand side` measured **by** `right-hand side`,
#'     }
#'     \item{on}{
#'       `left-hand side` regressed **on** `right-hand side`,
#'     }
#'     \item{with}{
#'       `left-hand side` covarying **with** `right-hand side`,
#'     }
#'     \item{on 1}{
#'       `left-hand side` regressed **on 1** for mean structure.
#'     }
#'   }
#'
#' @return \eqn{n} variates from \eqn{ \mathbf{X} \sim \mathcal{N}_{k}
#'     \left(
#'       \boldsymbol{\mu} \left( \boldsymbol{\theta} \right),
#'       \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     \right) .
#'   }
#'
#' @section par.label:
#'   Each parameter should have a numeric value.
#'
#' @section Line breaks:
#'   The characters `\n` and `;` can be used as line breaks.
#'   **Each line should end with a line break.**
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family data generation functions
#' @keywords data
#'
#' @inherit ramR references
#'
#' @inheritParams Eq2RAM
#' @inheritParams RAM2Data
#'
#' @examples
#' eq <- "
#'   # lhs op   rhs value
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' Eq2Data(eq, n = 100)
#' @export
Eq2Data <- function(eq,
                    n,
                    ...) {
  Expectations <- Eq2RAM(
    eq,
    par = FALSE
  )
  stopifnot(
    is.numeric(
      Expectations$par.table$par.label
    )
  )
  Expectations <- Expectations.default(
    A = Expectations$A,
    S = Expectations$S,
    u = Expectations$u,
    Filter = Expectations$Filter
  )
  return(
    MASS::mvrnorm(
      n = n,
      mu = Expectations$g,
      Sigma = Expectations$M,
      ...
    )
  )
}
