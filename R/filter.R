#' Generating group fixed-X knockoff variables
#'
#' This function samples multivariate Gaussian fixed-X knockoff variables
#'
#' @param x A \code{n}-by-\code{p} matrix of original variables
#' @param y The \code{n}-dimensional response vector. If the error standard error \code{sigma} is not specified, then \code{y} is required to estimate \code{sigma}.
#' @param group A \code{p}-dimensional vector of consecutive integers describing the group information of the variables. For example, group <- c(1, 1, 2, 2, 3, 3) means the first two vars belong to group 1, second two vars belong to group 2, and the last two vars belong to group 3. The 0-th group are fixed and not to be selected.
#' @param method Method of constructing fixed-X knockoffs.
#' @param loss a character string specifying the loss function to use in fitting the group lasso problem. Default is \code{"ls"} for least squares loss.
#' @param q A pre-specified nominal level for the controlled group FDR
#' @param plus Default to be \code{TRUE}, which gives the knockoff+ procedure that controls the FDR under \code{q}. If \code{plus = FALSE}, then the output does not technically controls FDR under \code{q}, but it controls the modified group FDR under \code{q}.
#' @param sigma The error standard deviation. If specified as \code{NULL} (as default), then the program will estimate it. In such a case, the response vector \code{y} has to be supplied to the program.
#' @param randomize An indicator of whether the knockoffs are constructed deterministically (default) or randomized.
#'
#' @return An object of class "gknockoff.result", which contains the following components:
#' \describe{
#' \item{\code{call}}{function call}
#' \item{\code{x}}{the \code{n}-by-\code{p} original design matrix. This could be transformed or augmented in the case where \code{n} < \code{2p}}
#' \item{\code{xk}}{the \code{n}-by-\code{p} constructed group fixed-X knockoff variables}
#' \item{\code{y}}{the \code{n}-dim response vector, if specified into the program. This could be the augmented if \code{n} < \code{2p}}
#' \item{\code{group}}{grouping information}
#' \item{\code{stat}}{group knockoff statistics}
#' \item{\code{threshold}}{the computed thresholed for filter}
#' \item{\code{selected}}{the set of selected groups}
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' p <- 10
#' gr <- c(rep(0, 2), rep(1, 3), rep(2, 2), 3, rep(4, 2))
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- rep(0, p)
#' beta[gr == 0] <- 10
#' beta[gr == 2] <- 2
#' beta[gr == 3] <- 5
#' y <- x %*% beta + rnorm(n)
#' selected <- filter_gko(x = x, y  = y, group = gr)
#' @export
#'
filter_gko <- function(x, y, group, method = "equi", loss = "ls", q = 0.2, plus = TRUE, sigma = NULL, randomize = FALSE){
  # we first grab the variables that need to be selected
  xm <- x[, which(group != 0)]
  gm <- group[group != 0]
  # offset: variables that are known to be in the model, and thus will not be selected
  xo <- x[, which(group == 0)]

  # construct group knockoffs
  # gkoo <- knockoff::create.fixed(X = xm, method = "equi")
  gko <- create_fix_gko(x = xm, y = y, group = gm, method = method, sigma = sigma, randomize = randomize)

  # compute group knockoff statistics
  stat <- stat_gglasso_lambda(xo = xo, gko = gko, group = group, loss = loss)

  # compute the filtering threshold
  t <- thresh_gko(stat = stat, q = q, plus = plus)
  # return the selected groups
  selected <- which(stat >= t)

  return(structure(list(call = match.call(), xo = xo, x = gko$x, xk = gko$xk, y = gko$y,
                        group = group, stat = stat, threshold = t, selected = selected),
                   class = "gknockoff.result"))
}
