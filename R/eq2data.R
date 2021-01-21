#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Equations to Sample Data
#'
#' @description Generates data from a multivariate normal distribution
#'   from an input model with specified values.
#'
#' @details The input model is a character string
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
#' @inherit ramR references
#' @inheritParams mvn
#' @param model Character string. Input model. See Details.
#' @examples
#' model <- "
#'   # VARIABLE1 OPERATION VARIABLE2 VALUE
#'   e           by        y         1.00;
#'   y           on        x         1.00;
#'   e           with      e         0.25;
#'   x           with      x         0.25;
#'   y           on        1         0.00;
#'   x           on        1         0.50
#' "
#' eq2data(model, n = 100)
#' @export
eq2data <- function(model,
                    n,
                    ...) {
  RAM <- eq2ram(model)
  if (is.null(RAM$u)) {
    stop(
      "Mean structure is required."
    )
  }
  if (!is.numeric(RAM$A)) {
    stop(
      "The `value` column in the input model should be numeric."
    )
  }
  out <- mvn(
    n = n,
    A = RAM$A,
    S = RAM$S,
    u = RAM$u,
    filter = RAM$filter,
    ...
  )
  model <- RAM$model
  colnames(model) <- c("var1", "op", "var2", "value")
  attributes(out)$model <- model
  return(
    out
  )
}
