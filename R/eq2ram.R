#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Equations to RAM Matrices
#'
#' @description Converts model equations to RAM matrices.
#'
#' @details The input model is a character string
#'   that specifies the associations between the variables.
#'
#' @section Syntax:
#'   Each line should follow the syntax below
#'
#'   `VARIABLE1 OPERATION VARIABLE2 LABEL`
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
#' @inherit ramR references
#' @param model Character string. Input model. See Details.
#' @return Returns the input model
#'   and the `A`, `S`, `filter`, and `u` RAM matrices.
#' @examples
#' model <- "
#'   # VARIABLE1 OPERATION VARIABLE2 LABEL
#'   e           by        y         1;
#'   y           on        x         beta;
#'   e           with      e         sigma[varepsilon]^2;
#'   x           with      x         sigma[x]^2;
#'   y           on        1         alpha;
#'   x           on        1         mu[x]
#' "
#' eq2ram(model)
#' @export
eq2ram <- function(model) {
  model <- gsub(pattern = "#[^\\\n]*", replacement = "", x = model)
  model <- gsub(pattern = "\n", replacement = "", model)
  model <- unlist(strsplit(x = model, split = ";"))
  model <- trimws(x = gsub(pattern = "\\s+", replacement = " ", x = model))
  model <- do.call(what = "rbind", args = strsplit(x = model, split = " "))
  colnames(model) <- c("var1", "op", "var2", "label")
  by <- model[which(model[, "op"] == "by"), , drop = FALSE]
  with <- model[which(model[, "op"] == "with"), , drop = FALSE]
  on <- model[which(model[, "op"] == "on" & model[, "var2"] != "1"), , drop = FALSE]
  one <- model[which(model[, "op"] == "on" & model[, "var2"] == "1"), , drop = FALSE]
  v <- unique(c(model[, "var1"], model[, "var2"]))
  v <- v[which(v != "1")]
  t <- length(v)
  h <- unique(by[, "var1"])
  g <- v[!v %in% h]
  v <- c(g, h)
  q <- length(h)
  p <- t - q
  A <- S <- matrix(
    data = 0,
    nrow = t,
    ncol = t
  )
  filter <- matrix(
    data = 0,
    nrow = p,
    ncol = t
  )
  diag(filter) <- 1
  colnames(filter) <- rownames(A) <- colnames(A) <- rownames(S) <- colnames(S) <- v
  rownames(filter) <- g
  # loadings
  for (i in seq_len(dim(by)[1])) {
    loadings <- by[i, , drop = FALSE]
    A[loadings[, "var2"], loadings[, "var1"]] <- to.numeric(loadings[, "label"])
  }
  # regressions
  for (i in seq_len(dim(on)[1])) {
    regressions <- on[i, , drop = FALSE]
    A[regressions[, "var1"], regressions[, "var2"]] <- to.numeric(regressions[, "label"])
  }
  # variances
  for (i in seq_len(dim(with)[1])) {
    variances <- with[i, , drop = FALSE]
    S[variances[, "var1"], variances[, "var2"]] <- to.numeric(variances[, "label"])
    S[variances[, "var2"], variances[, "var1"]] <- to.numeric(variances[, "label"])
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
      u[means[, "var1"], 1] <- to.numeric(means[, "label"])
    }
  } else {
    u <- NULL
  }
  model <- as.data.frame(model)
  label <- sapply(
    X = model$label,
    FUN = to.numeric
  )
  model$label <- label
  return(
    list(
      model = model,
      A = A,
      S = S,
      filter = filter,
      u = u
    )
  )
}

to.numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) {
    return(x)
  } else {
    return(out)
  }
}
