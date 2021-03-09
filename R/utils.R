#' AsNumeric
#'
#' By default, the function converts the input to numeric
#' if it can be coerced.
#' Otherwise, the original input is returned.
#' If all elements in the vector can be coerced to numeric,
#' the resulting vector will be numeric.
#'
#' @family utility functions
#' @keywords utils
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Vector.
#' @param return_NA Logical.
#'   If `return_NA = TRUE`,
#'   returns `NA` when `x` cannot be coerced.
#'   If `return_NA = FALSE`,
#'   returns the original item when `x` cannot be coerced.
#' @examples
#' x <- c("1", "2", "3")
#' AsNumeric(x)
#' x <- c("1", "2", "3", "a")
#' AsNumeric(x, return_NA = FALSE)
#' AsNumeric(x, return_NA = TRUE)
#' @export
AsNumeric <- function(x, return_NA = FALSE) {
  foo <- function(x, return_NA) {
    num <- suppressWarnings(
      as.numeric(x)
    )
    if (return_NA) {
      return(num)
    } else {
      if (is.na(num)) {
        return(x)
      } else {
        return(num)
      }
    }
  }
  return(
    sapply(
      X = x,
      FUN = foo,
      return_NA = return_NA
    )
  )
}

Length <- function(x,
                   transpose = FALSE) {
  if (transpose) {
    return(
      as.numeric(
        Ryacas::yac_str(
          paste0(
            "Length(",
            "Transpose(",
            x,
            ")",
            ")"
          )
        )
      )
    )
  } else {
    return(
      as.numeric(
        Ryacas::yac_str(
          paste0(
            "Length(",
            x,
            ")"
          )
        )
      )
    )
  }
}
