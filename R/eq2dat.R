#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Equations to Sample Data
#'
#' @description Generates data from a multivariate normal distribution
#'   from model equations.
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
#'   **Each line should end with a semicolon (`;`).**
#'
#' @section Value:
#'   Each parameter should have a numeric value.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @inherit eq2exp_num references
#' @inheritParams eq2exp_num
#' @inheritParams ram2dat
#' @examples
#' eq <- "
#'   # VARIABLE1 OPERATION VARIABLE2 VALUE
#'   e           by        y         1.00;
#'   y           on        x         1.00;
#'   e           with      e         0.25;
#'   x           with      x         0.25;
#'   y           on        1         0.00;
#'   x           on        1         0.50
#' "
#' eq2dat(eq, n = 100)
#' @export
eq2dat <- function(eq,
                   n,
                   ...) {
  ram <- eq2ram(eq)
  if (is.null(ram$u)) {
    stop(
      "Mean structure is required."
    )
  }
  return(
    ram2dat(
      n = n,
      A = ram$A,
      S = ram$S,
      u = ram$u,
      filter = ram$filter,
      ...
    )
  )
}
