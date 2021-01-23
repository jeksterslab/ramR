#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Equations to Expectations - Numeric
#'
#' @description Converts equations to expectations.
#'
#' @details The input is a character string
#'   that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `VARIABLE1 OPERATION VARIABLE2 VALUE`
#'
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{`VARIABLE1` measured **by** `VARIABLE2`}
#'     \item{on}{`VARIABLE1` regressed **on** `VARIABLE2`}
#'     \item{with}{`VARIABLE1` covarying **with** `VARIABLE2`}
#'     \item{on 1}{`VARIABLE1` regressed **on 1** for mean structure}
#'   }
#'
#' @section Value:
#'   Each parameter should have a numeric value.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @inherit eq2ram references
#' @inheritParams eq2ram
#' @examples
#' eq <- "
#'   # VARIABLE1 OPERATION VARIABLE2 VALUE
#'   e           by        y         1.00
#'   y           on        x         1.00
#'   e           with      e         0.25
#'   x           with      x         0.25
#'   y           on        1         0.00
#'   x           on        1         0.50
#' "
#' eq2exp_num(eq)
#' @export
eq2exp_num <- function(eq) {
  ram <- eq2ram(eq)
  variables <- list(
    variables = ram[["variables"]]
  )
  ram <- ram_num(
    A = ram$A,
    S = ram$S,
    u = ram$u,
    filter = ram$filter
  )
  return(
    c(
      variables,
      ram
    )
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Equations to Expectations - Symbolic
#'
#' @description Converts equations to expectations.
#'
#' @details The input is a character string
#'   that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `VARIABLE1 OPERATION VARIABLE2 LABEL/VALUE`
#'
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{`VARIABLE1` measured **by** `VARIABLE2`}
#'     \item{on}{`VARIABLE1` regressed **on** `VARIABLE2`}
#'     \item{with}{`VARIABLE1` covarying **with** `VARIABLE2`}
#'     \item{on 1}{`VARIABLE1` regressed **on 1** for mean structure}
#'   }
#'
#' @section Label:
#'   Each parameter should be labeled.
#'   The `LABEL` should be a number for fixed parameters
#'   and a variable name for free parameters.
#'   Equality contraints can be imposed by using the same variable name.
#'   Variable names can be written
#'   following the `R` mathematical notation `grDevices::plotmath()`.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @inherit eq2ram references
#' @inheritParams eq2ram
#' @examples
#' eq <- "
#'   # VARIABLE1 OPERATION VARIABLE2 LABEL/VALUE
#'   e           by        y         1
#'   y           on        x         beta
#'   e           with      e         sigma[varepsilon]^2
#'   x           with      x         sigma[x]^2
#'   y           on        1         alpha
#'   x           on        1         mu[x]
#' "
#' eq2exp_sym(eq)
#' @export
eq2exp_sym <- function(eq) {
  ram <- eq2ram(eq)
  variables <- list(
    variables = ram[["variables"]]
  )
  ram <- ram_sym(
    A = ram$A,
    S = ram$S,
    u = ram$u,
    filter = ram$filter
  )
  return(
    c(
      variables,
      ram
    )
  )
}
