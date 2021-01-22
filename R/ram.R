#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Expectations from the Reticular Action Model (RAM) Matrices - Numeric
#'
#' @inheritParams M_num
#' @inheritParams g_num
#' @examples
#' # This is a numerical example for the model
#' # y = alpha + beta * x + e
#' # y = 0 + 1 * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, 1, 1)
#' diag(S) <- c(0, 0.25, 1)
#' u <- c(0.00, 0.50, 0.00)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' ram_num(A, S, u, filter)
#' @export
ram_num <- function(A,
                    S,
                    u = NULL,
                    filter) {
  if (!is.null(u)) {
    u <- as.matrix(u)
    v <- v_num(
      A,
      u
    )
    g <- filter %*% v
    colnames(g) <- "g"
  } else {
    v <- NULL
    g <- NULL
  }
  C <- C_num(
    A,
    S
  )
  M <- filter %*% C %*% t(filter)
  return(
    list(
      A = A,
      S = S,
      u = u,
      filter = filter,
      v = v,
      g = g,
      C = C,
      M = M
    )
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Expectations from the Reticular Action Model (RAM) Matrices - Symbolic
#'
#' @inheritParams M_num
#' @inheritParams g_num
#' @export
#' @inheritParams M_sym
#' @inheritParams g_sym
#' @examples
#' # This is a symbolic example for the model
#' # y = alpha + beta * x + e
#'
#' A <- S <- matrix(
#'   data = 0,
#'   nrow = 3,
#'   ncol = 3
#' )
#' A[1, ] <- c(0, "beta", 1)
#' diag(S) <- c(0, "sigma2x", "sigma2e")
#' u <- c("alpha", "mux", 0)
#' filter <- diag(2)
#' filter <- cbind(filter, 0)
#' ram_sym(A, S, u, filter)
#' @export
ram_sym <- function(A,
                    S,
                    u = NULL,
                    filter) {
  filter <- Ryacas::ysym(filter)
  if (!is.null(u)) {
    u <- as.matrix(u)
    v <- v_sym(
      A,
      u
    )
    g <- filter * v
    u <- Ryacas::ysym(u)
  } else {
    v <- NULL
    g <- NULL
  }
  C <- C_sym(
    A,
    S
  )
  M <- filter * C * t(filter)
  return(
    list(
      A = Ryacas::ysym(A),
      S = Ryacas::ysym(S),
      u = u,
      filter = filter,
      v = v,
      g = g,
      C = C,
      M = M
    )
  )
}
