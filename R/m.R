#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Mean Structure Vector \eqn{\mathbf{m}}
#'
#' @description Derives the mean structure vector \eqn{\mathbf{m}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The mean structure vector \eqn{\mathbf{m}}
#'   as a function of Reticular Action Model (RAM) matrices is given by
#'
#'   \deqn{
#'     \mathbf{m}
#'     =
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)
#'     \mathbf{F}^{\mathsf{T}}
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{F}_{j \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'     is the \eqn{t \times 1} model-implied mean vector
#'   - \eqn{j} number of observed variables,
#'   - \eqn{k} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{j + k} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams Sigmatheta
#' @inherit Sigmatheta references
#' @param mutheta `t x 1` numeric vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)_{t \times 1}} .
#'   Model-implied mean vector.
#' @return Returns the mean structure vector \eqn{\mathbf{m}}
#'   derived from the `mutheta`, `A`, and `filter` matrices.
#' @export
m <- function(mutheta,
              A,
              filter) {
  # I - A
  IminusA <- diag(nrow(A)) - A
  # F^T * mu(theta)
  Ftmutheta <- crossprod(
    x = filter,
    y = mutheta
  )
  return(
    # F * (I - A) * F^T * mu(theta)
    filter %*% IminusA %*% Ftmutheta
  )
}
