#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Mean Structure Vector \eqn{\mathbf{M}}
#'
#' @description Derives the mean structure vector \eqn{\mathbf{M}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The mean structure vector \eqn{\mathbf{M}}
#'   as a function of Reticular Action Model (RAM) matrices is given by
#'
#'   \deqn{
#'     \mathbf{M}
#'     =
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)
#'     \mathbf{F}^{\mathsf{T}}
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'     is the \eqn{k \times 1} model-implied mean vector
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams Sigmatheta
#' @inherit Sigmatheta references
#' @param mutheta `k x 1` numeric vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)_{k \times 1}} .
#'   Model-implied mean vector.
#' @return Returns the mean structure vector \eqn{\mathbf{M}}
#'   derived from the \eqn{\mathbf{A}}, \eqn{\mathbf{F}},
#'   \eqn{ \mathbf{I}}, matrices and
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   vector.
#' @export
M <- function(mutheta,
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
