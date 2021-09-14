#' Compute group lasso statistics for group knockoffs
#'
#' This function computes the group knockoff statistics used for the filter step based on the maximum value of lambdas for which the group of variables have non-zero coefficients estimates.
#'
#' This function computes and returns a vector \code{stat}
#' \deqn{stat_j = max(\lambda_j, \tau_j) * sgn(\lambda_j - \tau_j),}
#' where \eqn{\lambda_j} and \eqn{\tau_j} are the largest values of lambdas at which the j-th group of variables, and the corresponding group knockoff variables, enter the model (i.e., coefficients become nonzero), respectively.
#'
#' @param xo The offset variable matrix
#' @param gko An object from group-knockoff construction.
#' @param group A \code{p}-dimensional vector of consecutive integers describing the group information of the variables. For example, group <- c(1, 1, 2, 2, 3, 3) means the first two vars belong to group 1, second two vars belong to group 2, and the last two vars belong to group 3. The 0-th group are fixed and not to be selected.
#' @param loss a character string specifying the loss function to use in fitting the group lasso problem. Default is \code{"ls"} for least squares loss.
#'
#' @import gglasso
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
#' stat <- stat_gglasso_lambda(gko = gko, group = gr)
#' @export
#'
stat_gglasso_lambda <- function(xo = NULL, gko, group, loss = "ls"){
  # input check
  if (!requireNamespace("gglasso", quietly = T))
    stop("gglasso package is not installed", call. = F)

  d <- ncol(gko$x)

  # augmented design matrix
  xx <- as.matrix(cbind(xo, gko$x, gko$xk))

  # unique non-offset group
  gm <- group[group != 0]
  # unique (non-offset) groups
  ug <- sort(unique(gm))
  # group used in gglasso
  gr <- c(rep(0, sum(group == 0)), gm, gm + length(ug))
  # penalty factor
  pf <- sqrt(as.integer(as.numeric(table(gr))))

  # fit the group lasso with the augmented design matrix xx
  # note that group = 0 is not accepted by gglasso
  # so we increment all group number by 1
  fit <- gglasso::gglasso(x = xx, y = gko$y,
                          group = gr + 1, pf = pf, loss = loss)

  first_nonzero <- function(x) match(T, abs(x) > 0)
  # find, on the element-wise fashion, the first index of lambda list whose corresponding estimate is non-zero
  idx_all <- apply(fit$beta[which(gr != 0), ], 1, first_nonzero)
  names(idx_all) <- NULL

  # translate the element-wise index to group-wise index
  idx_o <- head(idx_all, n = d)
  idx_k <- tail(idx_all, n = d)

  gind_o <- rep(NA, length(ug))
  gind_k <- rep(NA, length(ug))
  for(i in ug){
    gind_o[i] <- min(idx_o[which(gm == i)])
    gind_k[i] <- min(idx_k[which(gm == i)])
  }

  # get the corresponding lambda values
  lam_o <- fit$lambda[gind_o]
  lam_o[is.na(lam_o)] <- 0
  lam_k <- fit$lambda[gind_k]
  lam_k[is.na(lam_k)] <- 0

  out <- pmax(lam_o, lam_k) * sign(lam_o - lam_k)
  return(out)
}
