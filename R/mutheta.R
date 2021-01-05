#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Model-Implied Mean Vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied mean vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied mean vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \mathbf{F}^{\mathsf{T}}
#'     \mathbf{M}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{M}_{m \times 1}} represents the mean structure,
#'     that is, a vector of means and intercepts,
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams Sigmatheta
#' @inherit Sigmatheta references
#' @param M `m x 1` numeric vector \eqn{\mathbf{M}_{m \times 1}}.
#'   Mean structure. Vector of means and intercepts.
#' @return Returns the model-implied mean vector
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   derived from the `M`, `A`, and `filter` matrices.
#' @export
mutheta <- function(M,
                    A,
                    filter) {
  # (I - A)^{-1}
  invIminusA <- solve(
    diag(nrow(A)) - A
  )
  # F^{T} * M
  FtM <- crossprod(
    x = filter,
    y = M
  )
  return(
    # F * (I - A)^{-1} * F^{T} * M
    filter %*% invIminusA %*% FtM
  )
}
