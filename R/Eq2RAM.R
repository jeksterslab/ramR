#' Equations to RAM Matrices
#'
#' Converts model equations to RAM matrices.
#'
#' The argument `eq` is a character string
#' that specifies the associations between the variables.
#' See `Syntax`, `Operations`, `par.label`,
#' `Line breaks`, and `Comments` below.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs <space> op <space> rhs <space> par.label <\n> or <;>`
#'
#'   \describe{
#'     \item{lhs}{
#'       is the variable on the **left-hand side**,
#'     }
#'     \item{rhs}{
#'       is the variable on the **right-hand side**,
#'     }
#'     \item{op}{
#'       is the **operation** between `lhs` and `rhs`,
#'     }
#'     \item{par.label}{
#'       is the column of **parameter label**,
#'     }
#'     \item{\\n or ;}{
#'       are **line breaks**.
#'       **Each line should end with a line break.**
#'     }
#'   }
#'
#' @section Operations:
#'   The associations are defined by the following operations
#'
#'   \describe{
#'     \item{by}{
#'       `left-hand side` measured **by** `right-hand side`,
#'     }
#'     \item{on}{
#'       `left-hand side` regressed **on** `right-hand side`,
#'     }
#'     \item{with}{
#'       `left-hand side` covarying **with** `right-hand side`,
#'     }
#'     \item{on 1}{
#'       `left-hand side` regressed **on 1** for mean structure.
#'     }
#'   }
#'
#' @section par.label:
#'   Each parameter should be labeled.
#'   The `par.label` should be a number for fixed parameters
#'   and a character string for free parameters.
#'   Equality contraints can be imposed by using the same `par.label`.
#'
#' @section Line breaks:
#'   The characters `\n` and `;` can be used as line breaks.
#'   **Each line should end with a line break.**
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @return Returns a list with the following elements
#'
#'   \describe{
#'     \item{par.table}{Parameter table.}
#'     \item{variables}{Variable names.}
#'     \item{g.variables}{Variable names of observed variables.}
#'     \item{h.variables}{Variable names of latent variables.}
#'     \item{A}{
#'       `t by t` matrix \eqn{\mathbf{A}}.
#'       Asymmetric paths (single-headed arrows),
#'       such as regression coefficients and factor loadings.
#'     }
#'     \item{S}{
#'       `t by t` numeric matrix \eqn{\mathbf{S}}.
#'       Symmetric paths (double-headed arrows),
#'       such as variances and covariances.
#'     }
#'     \item{u}{
#'       `t by 1` matrix \eqn{\mathbf{u}} of mean structure parameters.
#'     }
#'     \item{Filter}{
#'       `p by t` numeric matrix
#'       \eqn{\mathbf{F}}.
#'       Filter matrix used to select observed variables.
#'     }
#'   }
#'
#' `par.tables` in the list is a data.frame
#' with the following columns
#'
#'   \describe{
#'     \item{lhs}{
#'       is the variable on the **left-hand side**,
#'     }
#'     \item{rhs}{
#'       is the variable on the **right-hand side**,
#'     }
#'     \item{op}{
#'       is the **operation** between `lhs` and `rhs`,
#'     }
#'     \item{par.label}{
#'       is the column of **parameter label**,
#'     }
#'     \item{par.start}{
#'       is the column of **starting values** for estimation
#'       if `eq` has a fifth colulmn, and
#'     }
#'     \item{par.names}{
#'       is the column of **parameter label
#'       with `NAs` on fixed parameters**,
#'     }
#'     \item{par.type}{
#'       is the type of the parameter,
#'     }
#'     \item{RAM}{
#'       is the RAM matrix used to represent the parameter,
#'     }
#'     \item{RAM.row}{
#'       is the row index in the RAM matrix for the parameter, and
#'     }
#'     \item{RAM.col}{
#'       is the column index in the RAM matrix for the parameter.
#'     }
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family eq functions
#' @keywords eq
#'
#' @inherit ramR references
#'
#' @inheritParams EqParse
#' @param par Logical.
#'   If `par = TRUE`, use `par.index` as labels.
#'   If `par = FALSE`, use `par.label` as labels.
#'
#' @examples
#' # Numeric ----------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs par.label
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' Eq2RAM(eq)
#'
#' # Symbolic ----------------------------------------------------------
#' eq <- "
#'   # lhs op   rhs par.label
#'     e   by   y   1
#'     y   on   x   beta
#'     e   with e   sigmae2
#'     x   with x   sigmax2
#'     y   on   1   alpha
#'     x   on   1   mux
#' "
#' Eq2RAM(eq)
#' @export
Eq2RAM <- function(eq,
                   par = FALSE) {
  if (methods::is(eq, "ParameterTable")) {
    par.table <- eq
  } else {
    par.table <- EqParse(eq)
  }
  temp <- par.table
  if (par) {
    temp[, "par.label"] <- par.table[, "par.index"]
  }
  by <- temp[which(
    temp[, "op"] == "by"
  ), , drop = FALSE]
  with <- temp[which(
    temp[, "op"] == "with"
  ), , drop = FALSE]
  on <- temp[which(
    temp[, "op"] == "on" & temp[, "rhs"] != "1"
  ), , drop = FALSE]
  one <- temp[which(
    temp[, "op"] == "on" & temp[, "rhs"] == "1"
  ), , drop = FALSE]
  v <- unique(
    c(
      temp[, "lhs"],
      temp[, "rhs"]
    )
  )
  v <- v[which(v != "1")]
  t <- length(v)
  h <- unique(by[, "lhs"])
  g <- v[!(v %in% h)]
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
  colnames(Filter) <- v
  rownames(A) <- v
  colnames(A) <- v
  rownames(S) <- v
  colnames(S) <- v
  rownames(Filter) <- g
  # loadings
  for (i in seq_len(dim(by)[1])) {
    loadings <- by[i, , drop = FALSE]
    A[loadings[, "rhs"], loadings[, "lhs"]] <- to.numeric(
      loadings[, "par.label"]
    )
    by$par.type[i] <- "loading"
    by$RAM[i] <- "A"
    by$RAM.row[i] <- loadings[, "rhs"]
    by$RAM.col[i] <- loadings[, "lhs"]
  }
  # regressions
  for (i in seq_len(dim(on)[1])) {
    regressions <- on[i, , drop = FALSE]
    A[regressions[, "lhs"], regressions[, "rhs"]] <- to.numeric(
      regressions[, "par.label"]
    )
    on$par.type[i] <- "regression"
    on$RAM[i] <- "A"
    on$RAM.row[i] <- regressions[, "lhs"]
    on$RAM.col[i] <- regressions[, "rhs"]
  }
  # variances
  for (i in seq_len(dim(with)[1])) {
    variances <- with[i, , drop = FALSE]
    S[variances[, "lhs"], variances[, "rhs"]] <- to.numeric(
      variances[, "par.label"]
    )
    S[variances[, "rhs"], variances[, "lhs"]] <- to.numeric(
      variances[, "par.label"]
    )
    with$par.type[i] <- "variance/covariance"
    with$RAM[i] <- "S"
    with$RAM.row[i] <- variances[, "lhs"]
    with$RAM.col[i] <- variances[, "rhs"]
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
      u[means[, "lhs"], 1] <- to.numeric(
        means[, "par.label"]
      )
      one$par.type[i] <- "mean/intercept"
      one$RAM[i] <- "u"
      one$RAM.row[i] <- means[, "lhs"]
      one$RAM.col[i] <- 1
    }
  } else {
    u <- NULL
  }
  # RAM matrices indices
  temp <- rbind(
    by,
    on,
    with,
    one
  )
  for (i in seq_len(nrow(temp))) {
    for (j in seq_along(v)) {
      if (temp$RAM.row[i] == v[j]) {
        temp$RAM.row[i] <- j
      }
      if (temp$RAM.col[i] == v[j]) {
        temp$RAM.col[i] <- j
      }
    }
  }
  temp <- temp[order(as.numeric(rownames(temp))), ]
  par.table$par.type <- temp$par.type
  par.table$RAM <- temp$RAM
  par.table$RAM.row <- temp$RAM.row
  par.table$RAM.col <- temp$RAM.col
  class(par.table) <- c("ParameterTable", class(par.table))
  out <- list(
    par.table = par.table,
    variables = v,
    g.variables = g,
    h.variables = h,
    A = A,
    S = S,
    u = u,
    Filter = Filter
  )
  class(out) <- c("RAM", class(out))
  return(out)
}
