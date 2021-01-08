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
#'     \mathbf{F}
#'     \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \boldsymbol{\Omega}
#'     \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{\mathsf{T}}
#'     \mathbf{F}^{\mathsf{T}}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{t \times t}} represents asymmetric paths
#'     (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\boldsymbol{\Omega}_{t \times t}} represents symmetric paths
#'     (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{F}_{j \times t}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{I}_{t \times t}} represents an identity matrix,
#'   - \eqn{j} number of observed variables,
#'   - \eqn{k} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{j + k} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A `t x t` numeric matrix
#'   \eqn{\mathbf{A}_{t \times t}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param Omega `t x t` numeric matrix
#'   \eqn{\boldsymbol{\Omega}_{t \times t}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @param filter `j x t` numeric matrix
#'   \eqn{\mathbf{F}_{j \times t}}.
#'   Filter matrix used to select variables.
#'   If `filter = NULL`,
#'   the filter matrix \eqn{\mathbf{F}} is assumed to be
#'   a \eqn{t \times t} identity matrix.
#' @return Returns the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   derived from the `A`, `Omega`, and `filter` matrices.
#' @references
#'   McArdle, J. J. (2005).
#'   The development of the RAM rules
#'   for latent variable structural equation modeling.
#'   In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#'   *Contemporary Psychometrics: A festschrift for Roderick P. McDonald*
#'   (pp. 225--273).
#'   Lawrence Erlbaum Associates.
#'
#'   McArdle, J. J., & McDonald, R. P. (1984).
#'   Some algebraic properties
#'   of the Reticular Action Model for moment structures.
#'   *British Journal of Mathematical and Statistical Psychology*, *37* (2),
#'   234--251.
#' @export
Sigmatheta <- function(A,
                       Omega,
                       filter = NULL) {
  # (I - A)^{-1}
  invIminusA <- solve(
    diag(nrow(A)) - A
  )
  # Omega * ((I - A)^{-1})^T
  OmegaTinvIminusA <- tcrossprod(
    x = Omega,
    y = invIminusA
  )
  # (I - A)^{-1} * Omega * ((I - A)^{-1})^T
  inner <- invIminusA %*% OmegaTinvIminusA
  if (is.null(filter)) {
    return(inner)
  } else {
    return(
      # F * (I - A)^{-1} * Omega * ((I - A)^{-1})^T * F^T
      # filter %*% inner %*% t(filter)
      filter %*% tcrossprod(
        x = inner,
        y = filter
      )
    )
  }
}
