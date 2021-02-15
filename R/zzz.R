to.numeric <- function(x) {
  out <- suppressWarnings(
    as.numeric(x)
  )
  if (is.na(out)) {
    return(x)
  } else {
    return(out)
  }
}

isFALSE <- function(x) {
  return(
    identical(
      x,
      FALSE
    )
  )
}

R2Yac <- function(x) {
  if (methods::is(x, "yac_symbol")) {
    return(x)
  } else {
    return(Ryacas::ysym(x))
  }
}

RMatrix2Yac <- function(x) {
  x <- R2Yac(x)
  y_res <- Ryacas::yac_str(x$yacas_cmd)
  y <- Ryacas::ysym(y_res)
  stopifnot(y$is_mat)
  return(x)
}

RVector2Yac <- function(x,
                        col = FALSE) {
  x <- R2Yac(x)
  y_res <- Ryacas::yac_str(x$yacas_cmd)
  y <- Ryacas::ysym(y_res)
  if (col) {
    if (y$is_vec) {
      z <- x$yacas_cmd
      z <- gsub("\\{", "", z)
      z <- gsub("\\}", "", z)
      z <- unlist(strsplit(z, ","))
      z <- as.matrix(z)
      return(Ryacas::ysym(z))
    }
  }
  return(x)
}
