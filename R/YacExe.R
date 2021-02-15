#' Execute Yacas Expression using Ryacas
#'
#' Execute `yacas` expressions using `Ryacas`.
#'
#' @references
#' [YACAS - Yet Another Computer Algebra System](http://www.yacas.org/)
#'
#' [Ryacas](http://r-cas.github.io/ryacas/)
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family Yacas functions
#' @keywords yac
#'
#' @param expr Character string. Yacas expresion.
#' @param str Logical.
#'   If `str = TRUE`, returns symbolic result as string.
#'   If `str = FALSE`, returns symbolic result as an `R` expression.
#' @param ysym Logical.
#'   If `ysym = TRUE`, returns symbolic result as `yac_symbol` if `str = TRUE`.
#'   If `ysym = FALSE`, returns symbolic result as string if `str = TRUE`.
#' @param tex Logical.
#'   If `tex = TRUE` returns results as latex math if `str = TRUE`.
#'   If `tex = FALSE`, returns symbolic result as string if `str = TRUE`.
#' @param simplify Logical. Simplify symbolic results.
#' @examples
#' A <- Ryacas::ysym(
#'   matrix(c("a", "c", "b", "d"), ncol = 2)
#' )
#' expr <- paste0("Determinant(", A, ")")
#' YacExe(expr)
#' YacExe(expr, str = TRUE, ysym = FALSE)
#' YacExe(expr, str = TRUE, tex = TRUE)
#' YacExe(expr, str = FALSE)
#'
#' a <- 1
#' b <- 2
#' c <- 3
#' d <- 4
#' eval(YacExe(expr, str = FALSE))
#'
#' A <- matrix(c(a, c, b, d), ncol = 2)
#' det(A)
#' @export
YacExe <- function(expr,
                   str = TRUE,
                   ysym = TRUE,
                   tex = FALSE,
                   simplify = FALSE) {
  if (simplify) {
    expr <- paste0(
      "Simplify(", expr, ")"
    )
  }
  if (str) {
    out <- Ryacas::yac_str(expr)
    if (tex) {
      return(
        Ryacas::tex(
          Ryacas::ysym(out)
        )
      )
    } else {
      if (ysym) {
        return(
          Ryacas::ysym(out)
        )
      } else {
        return(
          out
        )
      }
    }
  } else {
    return(
      Ryacas::yac_expr(expr)
    )
  }
}
