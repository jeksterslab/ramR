#' Equations to RAM Matrices
#'
#' Converts model equations to RAM matrices.
#'
#' The input is a character string
#'   that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs operation rhs label`
#'
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{`left-hand side` measured **by** `right-hand side`}
#'     \item{on}{`left-hand side` regressed **on** `right-hand side`}
#'     \item{with}{`left-hand side` covarying **with** `right-hand side`}
#'     \item{on 1}{`left-hand side` regressed **on 1** for mean structure}
#'   }
#'
#' @section label:
#'   Each parameter should be labeled.
#'   The `label` should be a number for fixed parameters
#'   and a character string for free parameters.
#'   Equality contraints can be imposed by using the same label.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inherit ramR references
#' @param eq Character string. Equations. See Details.
#' @param par Logical.
#'   If `par = TRUE`, use `par.index` as labels.
#'   If `par = FALSE`, use `label` as labels.
#' @examples
#' # Numeric ----------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs label
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' Eq2RAM(eq, par = FALSE)
#' Eq2RAM(eq, par = TRUE)
#'
#' # Symbolic ----------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs label
#'     e   by   y   1
#'     y   on   x   beta
#'     e   with e   sigmae2
#'     x   with x   sigmax2
#'     y   on   1   alpha
#'     x   on   1   mux
#' "
#' Eq2RAM(eq, par = FALSE)
#' Eq2RAM(eq, par = TRUE)
#' @export
Eq2RAM <- function(eq,
                   par = FALSE) {
  eq <- Parse(eq)
  if (par) {
    eq[, "label"] <- eq[, "par.index"]
  }
  by <- eq[which(eq[, "op"] == "by"), , drop = FALSE]
  with <- eq[which(eq[, "op"] == "with"), , drop = FALSE]
  on <- eq[which(eq[, "op"] == "on" & eq[, "rhs"] != "1"), , drop = FALSE]
  one <- eq[which(eq[, "op"] == "on" & eq[, "rhs"] == "1"), , drop = FALSE]
  v <- unique(c(eq[, "lhs"], eq[, "rhs"]))
  v <- v[which(v != "1")]
  t <- length(v)
  h <- unique(by[, "lhs"])
  g <- v[!v %in% h]
  v <- c(g, h)
  q <- length(h)
  p <- t - q
  A <- S <- matrix(
    data = 0,
    nrow = t,
    ncol = t
  )
  Filter <- matrix(
    data = 0,
    nrow = p,
    ncol = t
  )
  diag(Filter) <- 1
  colnames(Filter) <- rownames(A) <- colnames(A) <- rownames(S) <- colnames(S) <- v
  rownames(Filter) <- g
  # loadings
  for (i in seq_len(dim(by)[1])) {
    loadings <- by[i, , drop = FALSE]
    A[loadings[, "rhs"], loadings[, "lhs"]] <- to.numeric(loadings[, "label"])
  }
  # regressions
  for (i in seq_len(dim(on)[1])) {
    regressions <- on[i, , drop = FALSE]
    A[regressions[, "lhs"], regressions[, "rhs"]] <- to.numeric(regressions[, "label"])
  }
  # variances
  for (i in seq_len(dim(with)[1])) {
    variances <- with[i, , drop = FALSE]
    S[variances[, "lhs"], variances[, "rhs"]] <- to.numeric(variances[, "label"])
    S[variances[, "rhs"], variances[, "lhs"]] <- to.numeric(variances[, "label"])
  }
  # means
  if (dim(one)[1] > 0) {
    u <- matrix(
      data = 0,
      nrow = t,
      ncol = 1
    )
    rownames(u) <- v
    colnames(u) <- "u"
    for (i in seq_len(dim(one)[1])) {
      means <- one[i, , drop = FALSE]
      u[means[, "lhs"], 1] <- to.numeric(means[, "label"])
    }
  } else {
    u <- NULL
  }
  # -------------------------------------------------------------
  return(
    list(
      eq = eq,
      variables = v,
      A = A,
      S = S,
      Filter = Filter,
      u = u
    )
  )
}
