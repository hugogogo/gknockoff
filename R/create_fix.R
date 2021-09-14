#' Generating group fixed-X knockoff variables
#'
#' This function samples multivariate Gaussian fixed-X knockoff variables
#'
#' @param x A \code{n}-by-\code{p} matrix of original variables
#' @param y The \code{n}-dimensional response vector. If the error standard error \code{sigma} is not specified, then \code{y} is required to estimate \code{sigma}.
#' @param group A \code{p}-dimensional vector of consecutive integers describing the group information of the variables. For example, group <- c(1, 1, 2, 2, 3, 3) means the first two vars belong to group 1, second two vars belong to group 2, and the last two vars belong to group 3
#' @param method Method of constructing fixed-X knockoffs.
#' @param sigma The error standard deviation. If specified as \code{NULL} (as default), then the program will estimate it. In such a case, the response vector \code{y} has to be supplied to the program.
#' @param randomize An indicator of whether the knockoffs are constructed deterministically (default) or randomized.
#'
#' @return An object of class "gknockoff.variables", which contains the following components:
#' \describe{
#' \item{\code{x}}{the \code{n}-by-\code{p} original design matrix. This could be transformed or augmented in the case where \code{n} < \code{2p}}
#' \item{\code{xk}}{the \code{n}-by-\code{p} constructed group fixed-X knockoff variables}
#' \item{\code{y}}{the \code{n}-dim response vector, if specified into the program. This could be the augmented if \code{n} < \code{2p}.}
#' }
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
#' @export
#'
create_fix_gko <- function(x, y = NULL, group, method = c("equi", "sdp"),
                           sigma = NULL, randomize = FALSE){
  method = match.arg(method)

  # xo are offset variables
  xo <- x[, which(group == 0)]
  # xm are variables to be selected
  xm <- x[, which(group != 0)]

  n <- nrow(xm)
  p <- ncol(xm)
  q <- ncol(xo)

  # only study the non-conditioning groups
  grm <- group[group != 0]
  # re-order the columns of xm
  # such that the columns are in the same group with increasing group number
  xm <- xm[, order(grm)]
  # xm[, order(order(grm))] will just give the very original xm back

  # gr gives grouping information in the ordered columns
  gr <- sort(grm)

  # now we look at the some pre-processing of x
  # i.e., if p < n < 2p, then we need to augment both x (with extra rows) and y
  if (n <= p)
    stop("Input non-conditioning design matrix must have dimensions n > p")
  if (n < 2 * p) {
    warning("Input non-conditioning design matrix has dimensions p < n < 2p. ", "Augmenting the model with extra rows.", immediate. = T)
    xm.add <- matrix(0, 2 * p - n, p)
    colnames(xm.add) <- colnames(xm)
    xm = rbind(xm, xm.add)

    if (is.null(sigma)) {
      if (is.null(y)) {
        stop("Either the noise level \"sigma\" or the response variables \"y\" must\n
             be provided in order to augment the data with extra rows.")
      }
      else {
        x.svd = svd(x, nu = n, nv = 0)
        u2 = x.svd$u[, (p + 1):n]
        sigma = sqrt(mean((t(u2) %*% y)^2))
      }
    }
    if (randomize)
      y.extra = rnorm(2 * p - n, sd = sigma)
    else y.extra = with_seed(0, rnorm(2 * p - n, sd = sigma))
    y = c(y, y.extra)
  }

  # now n >= 2p
  # make sure x has column mean zero and column norm 1
  x = standardize(x, center = FALSE)
  # now construct group knockoffs
  xk = switch(method,
              equi = create_equivariant(x, gr, randomize),
              sdp = create_sdp(x, gr, randomize))

  if(!is.null(colnames(xk)))
    colnames(xk) <- paste0(colnames(x), "_ko")
  # give the original order in the columns of x and xk
  x <- x[, order(order(group))]
  xk <- xk[, order(order(group))]
  structure(list(x = x, xk = xk, y = y), class = "gknockoff.variables")
}
