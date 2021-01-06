#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Matrix of symmetric paths (double-headed arrows)
#'   \eqn{\boldsymbol{\Omega}}
#'
#' @description Derives the matrix of symmetric paths (double-headed arrows)
#'   \eqn{\boldsymbol{\Omega}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The matrix of symmetric paths (double-headed arrows)
#'   \eqn{\boldsymbol{\Omega}}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \boldsymbol{\Omega}
#'     =
#'     \left( \mathbf{I} - \mathbf{A} \right)
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     \left( \mathbf{I} - \mathbf{A} \right)^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\boldsymbol{\Sigma}_{t \times t}} represents
#'     the model-implied variance-covariance matrix,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{j} number of observed variables,
#'   - \eqn{k} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{j + k} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams Sigmatheta
#' @inherit Sigmatheta references
#' @param Sigmatheta `t x t` numeric matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}.
#'   Model-implied variance-covariance matrix.
#' @return Returns the matrix of symmetric paths (double-headed arrows)
#'   \eqn{\boldsymbol{\Omega}}
#'   derived from the `A` and `Sigmatheta` matrices.
#' @export
Omega <- function(A,
                  Sigmatheta) {
  # I - A
  IminusA <- diag(nrow(A)) - A
  # Sigmatheta * (I - A)^{T}
  SigmathetaIminusAt <- tcrossprod(
    x = Sigmatheta,
    y = IminusA
  )
  return(
    # (I - A) * Sigmatheta * (I - A)^{T}
    IminusA %*% SigmathetaIminusAt
  )
}
