#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Model-Implied Variance-Covariance Matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{S}
#'     \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{\mathsf{T}} \mathbf{F}^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{m \times m}} represents symmetric paths (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A `m x m` numeric matrix
#'   \eqn{\mathbf{A}_{m \times m}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param S `m x m` numeric matrix
#'   \eqn{\mathbf{S}_{m \times m}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @param filter `k x m` numeric matrix
#'   \eqn{\mathbf{F}_{k \times m}}.
#'   Filter matrix used to select variables.
#' @return Returns the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   derived from the `A`, `S`, and `filter` matrices.
#' @references
#'   McArdle, J. J. (2013).
#'   The development of the RAM rules for latent variable structural equation modeling.
#'   In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#'   *Contemporary Psychometrics: A festschrift for Roderick P. McDonald* (pp. 225--273).
#'   Lawrence Erlbaum Associates.
#'
#'   McArdle, J. J., & McDonald, R. P. (1984).
#'   Some algebraic properties of the Reticular Action Model for moment structures.
#'   *British Journal of Mathematical and Statistical Psychology*, *37* (2), 234--251.
#' @export
Sigmatheta <- function(A,
                       S,
                       filter) {
  # (I - A)^{-1}
  invIminusA <- solve(
    diag(nrow(A)) - A
  )
  # S * ((I - A)^{-1})^T
  STinvIminusA <- tcrossprod(
    x = S,
    y = invIminusA
  )
  # (I - A)^{-1} * S * ((I - A)^{-1})^T
  inner <- invIminusA %*% STinvIminusA
  return(
    # F * (I - A)^{-1} * S * ((I - A)^{-1})^T * F^T
    filter %*% crossprod(
      x = inner,
      y = filter
    )
  )
}
