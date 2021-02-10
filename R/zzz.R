to.numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) {
    return(x)
  } else {
    return(out)
  }
}

isFALSE <- function(x) {
  identical(x, FALSE)
}

#' @param expr Character string. Yacas expresion.
#' @param str Logical.
#' If `TRUE`, returns symbolic result as string.
#' If `FALSE`, returns symbolic result as an `R` expression.
#' @param ysym Logical.
#' If `TRUE`, returns symbolic result as `yac_symbol` if `str = TRUE`.
#' If `FALSE`, returns symbolic result as string if `str = TRUE`.
#' @param simplify Logical. Simplify symbolic results.
#' @param tex Logical. Return results as latex if `str = TRUE`.
.exe <- function(expr,
                 str = TRUE,
                 ysym = TRUE,
                 simplify = FALSE,
                 tex = FALSE) {
  if (simplify) {
    expr <- paste0("Simplify(", expr, ")")
  }
  if (str) {
    out <- Ryacas::yac_str(expr)
    if (tex) {
      return(
        Ryacas::tex(Ryacas::ysym(out))
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
    return(Ryacas::yac_expr(expr))
  }
}

Parse <- function(eq) {
  eq <- gsub(pattern = "#[^\\\n]*", replacement = "", x = eq)
  eq <- unlist(strsplit(x = eq, split = "[\n;]"))
  # eq <- unlist(strsplit(x = eq, split = ";"))
  eq <- trimws(x = gsub(pattern = "\\s+", replacement = " ", x = eq))
  eq <- do.call(what = "rbind", args = strsplit(x = eq, split = " "))
  if (dim(eq)[2] == 5) {
    colnames(eq) <- c("lhs", "op", "rhs", "label", "start")
  } else {
    colnames(eq) <- c("lhs", "op", "rhs", "label")
  }
  eq[, "op"] <- tolower(eq[, "op"])
  # par.index-----------------------------------------------------------------
  label <- as.vector(eq[, "label"])
  uniquelabel <- rep(x = NA, length = length(eq[, "label"]))
  for (i in seq_along(label)) {
    if (is.na(suppressWarnings(as.numeric(label[i])))) {
      uniquelabel[i] <- label[i]
    } else {
      uniquelabel[i] <- NA
    }
  }
  uniquelabel <- unique(uniquelabel[stats::complete.cases(uniquelabel)])
  index <- paste0("p", seq_along(uniquelabel))
  par.index <- rep(x = NA, length = length(label))
  for (i in seq_along(label)) {
    for (j in seq_along(uniquelabel)) {
      if (uniquelabel[j] == label[i]) {
        par.index[i] <- index[j]
      }
    }
  }
  for (i in seq_along(par.index)) {
    if (is.na(par.index[i])) {
      par.index[i] <- label[i]
    }
  }
  # if all items are numeric the columns will be numeric
  eq <- cbind(
    eq,
    par.index
  )
  eq <- as.data.frame(eq)
  label <- sapply(
    X = eq$label,
    FUN = to.numeric
  )
  eq$label <- label
  par.index <- sapply(
    X = eq$par.index,
    FUN = to.numeric
  )
  eq$par.index <- par.index
  if ("start" %in% colnames(eq)) {
    start <- sapply(
      X = eq$start,
      FUN = function(x) suppressWarnings(as.numeric(x))
    )
    eq$start <- start
  }
  # ---------------------------------------------------------------------------
  return(eq)
}
