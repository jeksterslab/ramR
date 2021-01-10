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
#'     \mathbf{M}
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
#'   - \eqn{\mathbf{M}_{t \times 1}} represents the mean structure,
#'     that is, a vector of means and intercepts,
#'   - \eqn{j} number of observed variables,
#'   - \eqn{k} number of latent variables, and
#'   - \eqn{t} number of observed and latent variables, that is \eqn{j + k} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams Sigmatheta
#' @inherit Sigmatheta references
#' @param M `t x 1` numeric vector \eqn{\mathbf{M}_{t \times 1}}.
#'   Mean structure. Vector of means and intercepts.
#' @return Returns the model-implied mean vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   derived from the `M`, `A`, and `filter` matrices.
#' @export
mutheta <- function(M,
                    A,
                    filter = NULL) {
  if (is.vector(M)) {
    rowlabels <- names(M)
    M <- matrix(
      data = M,
      ncol = 1
    )
    rownames(M) <- rowlabels
  }
  if (is.null(filter)) {
    filter <- diag(nrow(A))
    colnames(filter) <- colnames(A)
  }
  E <- E(
    A = A
  )
  return(
    filter %*% E %*% M
  )
}
