#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title \eqn{\mathbf{A}} Matrix from Linear Regression Model Parameters
#'
#' @description Derives the \eqn{\mathbf{A}} matrix from linear regression model parameters.
#'
#' @details The linear regression model is given by
#'   \deqn{
#'     \mathbf{y}
#'     =
#'     \mathbf{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{y}} is an \eqn{n \times 1} vector of observations
#'     on the regressand variable
#'   - the data matrix \eqn{\mathbf{X}}
#'     (also known as design matrix, model matrix or regressor matrix)
#'     is an \eqn{n \times k} matrix of \eqn{n} observations of \eqn{k} regressors,
#'     which includes a regressor whose value is 1 for each observation on the first column
#'   - \eqn{\boldsymbol{\beta}} is a \eqn{k \times 1} vector regression coefficients
#'   - \eqn{\boldsymbol{\varepsilon}} is the stochastic error term distributed around 0
#'     with constant variance \eqn{\left( \sigma^2 \right)} across values of \eqn{\mathbf{X}}
#'
#'   Parameters of the linear regression model include
#'
#'   - \eqn{\boldsymbol{\beta}} and
#'   - \eqn{\sigma^2}
#'
#' @param beta Numeric vector or \eqn{k \times 1} matrix.
#'   Regression coefficients.
#' @export
A_linreg <- function(beta) {
  beta <- as.vector(beta)
  k <- length(beta)
  A <- matrix(
    data = 0,
    nrow = k,
    ncol = k
  )
  A[1, 2:k] <- beta[-1]
  return(A)
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title \eqn{\boldsymbol{\Omega}} Matrix from Linear Regression Model Parameters
#'
#' @description Derives the \eqn{\boldsymbol{Omega}} matrix from linear regression model parameters.
#'
#' @inherit A_linreg details
#' @param sigma2 Numeric.
#'   \eqn{\sigma^2} is the variance of the stochastic error term
#'   \eqn{\boldsymbol{\varepsilon}}.
#' @param SigmaX Numeric matrix.
#'   Covariances of \eqn{\mathbf{X}_{n \times 2, \cdots, k}}.
#' @export
Omega_linreg <- function(sigma2,
                         SigmaX) {
  if (length(as.vector(SigmaX)) == 1) {
    SigmaX <- matrix(
      data = SigmaX,
      ncol = 1
    )
  }
  if (!matrixR::is_sym(SigmaX)) {
    stop(
      "`SigmaX` should be a symmetric matrix."
    )
  }
  Omega <- cbind(
    0,
    SigmaX
  )
  Omega <- rbind(
    0,
    Omega
  )
  Omega[1, 1] <- sigma2
  return(Omega)
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model (RAM) Matrices from Linear Regression Model Parameters
#'
#' @description Derives the Reticular Action Model (RAM) matrices from linear regression model parameters.
#'
#' @inherit A_linreg details
#' @inheritParams A_linreg
#' @inheritParams Omega_linreg
#' @param muX Numeric vector of \eqn{k \times 1} matrix.
#'   Column means of \eqn{\mathbf{X}_{n \times 2, \cdots, k}}.
#' @export
ram_linreg <- function(beta,
                       sigma2,
                       SigmaX,
                       muX = NULL) {
  if (length(as.vector(SigmaX)) == 1) {
    SigmaX <- matrix(
      data = SigmaX,
      ncol = 1
    )
  }
  p <- dim(SigmaX)[1]
  k <- p + 1
  if (!(k == length(as.vector(beta)))) {
    stop(
      "`beta` and `SigmaX` have incompatible dimensions."
    )
  }
  A <- A_linreg(
    beta = beta
  )
  Omega <- Omega_linreg(
    sigma2 = sigma2,
    SigmaX = SigmaX
  )
  if (!is.null(muX)) {
    M <- matrix(
      data = NA_real_,
      ncol = 1,
      nrow = k
    )
    muX <- as.vector(muX)
    if (!(p == length(muX))) {
      stop(
        "`muX` and `SigmaX` have incompatible dimensions."
      )
    }
    M[1, 1] <- as.vector(beta)[1]
    M[2:k, 1] <- muX
    mutheta <- mutheta(
      M = M,
      A = A
    )
  } else {
    M <- NULL
    mutheta <- NULL
  }
  filter <- diag(nrow(A))
  Sigmatheta <- Sigmatheta(
    A = A,
    Omega = Omega,
    filter = filter
  )
  return(
    list(
      A = A,
      Omega = Omega,
      M = M,
      filter = filter,
      mutheta = mutheta,
      Sigmatheta = Sigmatheta
    )
  )
}
