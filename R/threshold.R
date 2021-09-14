#' Compute filter threshold for group knockoffs
#'
#' This function computes the threshold of the group knockoff statistics used for the filter step
#'
#' @param stat A \code{m} dimensional vector of group knockoff statistics
#' @param q A pre-specified nominal level for the controlled group FDR
#' @param plus Default to be \code{TRUE}, which gives the knockoff+ procedure that controls the FDR under \code{q}. If \code{plus = FALSE}, then the output does not technically controls FDR under \code{q}, but it controls the modified group FDR under \code{q}.
#'
#' @return a threshold
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' p <- 10
#' gr <- c(rep(1, 5), rep(2, 2), 3, rep(4, 2))
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- rep(0, p)
#' beta[gr == 2] <- 2
#' beta[gr == 3] <- 5
#' y <- x %*% beta + rnorm(n)
#' gko <- create_fix_gko(x = x, y = y, group = gr)
#' stat <- stat_gglasso_lambda(x = gko$x, xk = gko$xk, y = gko$y, group = gr)
#' t <- thresh_gko(stat, q = 0.2)
#' @export
thresh_gko <- function(stat, q, plus = TRUE){
  # ts <- sort(c(0, unique(abs(stat[stat != 0]))))
  ts <- sort(unique(abs(stat[stat != 0])))
  ratio <- sapply(ts, function(t) (as.integer(plus) + sum(stat <= -t))/max(1, sum(stat >= t)))
  ok <- which(ratio <= q)
  return(ifelse(length(ok) > 0, ts[ok[1]], Inf))
}
