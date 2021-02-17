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
