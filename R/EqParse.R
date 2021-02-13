#' Equation Parser
#'
#' Parse equations and return a parameter table (`par.table`).
#'
#' The input is a character string
#' that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `lhs operation rhs par.label par.start`
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
#' @section par.label:
#'   Each parameter should be labeled.
#'   The `par.label` should be a number for fixed parameters
#'   and a character string for free parameters.
#'   Equality contraints can be imposed by using the same `par.label`.
#'
#' @section par.start:
#'   Numerical values as starting values for estimation.
#'   Starting values for fixed parameters should be `NA`.
#'   Starting values should be identical
#'   for parameters constrained to be equal,
#'   that is, parameters with the same `par.label`.
#'
#' @section Line breaks:
#'   `\n` and `;` can be used as line breaks.
#'
#' @section Comments:
#'   Comments can be written after a hash (`#`) sign.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @family eq functions
#' @keywords eq
#'
#' @param eq Character string. Equations. See Details.
#'
#' @examples
#' eq <- "
#'   # lhs op   rhs par.label
#'     e   by   y   1
#'     y   on   x   beta
#'     e   with e   sigmae2
#'     x   with x   sigmax2
#'     y   on   1   alpha
#'     x   on   1   mux
#' "
#' EqParse(eq)
#'
#' eq <- "
#'   # lhs op   rhs value
#'     e   by   y   1
#'     y   on   x   1
#'     e   with e   1
#'     x   with x   0.25
#'     y   on   1   0
#'     x   on   1   0.50
#' "
#' EqParse(eq)
#'
#' eq <- "
#'   # lhs op   rhs par.label par.start
#'     e   by   y   1         NA
#'     y   on   x   beta      0.00
#'     e   with e   sigmae2   1.25
#'     x   with x   sigmax2   0.25
#'     y   on   1   alpha     0.00
#'     x   on   1   mux       0.50
#' "
#' EqParse(eq)
#'
#' # \n as line breaks------------------------------------------------
#'
#' eq <- "
#'   e by y 1 NA \n y on x beta 0.00
#'   e with e sigmae2 1.25 \n x with x sigmax2 0.25
#'   y on 1 alpha 0.00 \n x on 1 mux 0.50
#' "
#' EqParse(eq)
#'
#' # ; as line breaks------------------------------------------------
#'
#' eq <- "
#'   e by y 1 NA; y on x beta 0.00
#'   e with e sigmae2 1.25; x with x sigmax2 0.25
#'   y on 1 alpha 0.00; x on 1 mux 0.50
#' "
#' EqParse(eq)
#' @export
EqParse <- function(eq) {
  par.table <- gsub(
    pattern = "#[^\\\n]*",
    replacement = "",
    x = eq
  )
  par.table <- unlist(
    strsplit(
      x = par.table,
      split = "[\n;]"
    )
  )
  par.table <- trimws(
    x = gsub(
      pattern = "\\s+",
      replacement = " ",
      x = par.table
    )
  )
  par.table <- do.call(
    what = "rbind",
    args = strsplit(
      x = par.table,
      split = " "
    )
  )
  if (dim(par.table)[2] == 5) {
    colnames(par.table) <- c(
      "lhs",
      "op",
      "rhs",
      "par.label",
      "par.start"
    )
  } else {
    colnames(par.table) <- c(
      "lhs",
      "op",
      "rhs",
      "par.label"
    )
  }
  par.table[, "op"] <- tolower(
    par.table[, "op"]
  )
  # par.index-----------------------------------------------------------
  par.label <- as.vector(
    par.table[, "par.label"]
  )
  unique.par.label <- rep(
    x = NA,
    length = length(
      par.table[, "par.label"]
    )
  )
  for (i in seq_along(par.label)) {
    if (
      is.na(
        suppressWarnings(
          as.numeric(par.label[i])
        )
      )
    ) {
      unique.par.label[i] <- par.label[i]
    } else {
      unique.par.label[i] <- NA
    }
  }
  unique.par.label <- unique(
    unique.par.label[stats::complete.cases(unique.par.label)]
  )
  index <- paste0(
    "p",
    seq_along(unique.par.label)
  )
  par.index <- rep(
    x = NA,
    length = length(par.label)
  )
  for (i in seq_along(par.label)) {
    for (j in seq_along(unique.par.label)) {
      if (unique.par.label[j] == par.label[i]) {
        par.index[i] <- index[j]
      }
    }
  }
  for (i in seq_along(par.index)) {
    if (is.na(par.index[i])) {
      par.index[i] <- par.label[i]
    }
  }
  # if all items are numeric the columns will be numeric
  par.table <- cbind(
    par.table,
    par.index
  )
  par.table <- as.data.frame(
    par.table,
    stringsAsFactors = FALSE
  )
  par.label <- sapply(
    X = par.table$par.label,
    FUN = to.numeric
  )
  par.table$par.label <- par.label
  par.index <- sapply(
    X = par.table$par.index,
    FUN = to.numeric
  )
  par.table$par.index <- par.index
  if ("start" %in% colnames(par.table)) {
    start <- sapply(
      X = par.table$start,
      FUN = function(x) suppressWarnings(as.numeric(x))
    )
    par.table$start <- start
  }
  # --------------------------------------------------------------------
  return(par.table)
}
