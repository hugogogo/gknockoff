#' @export
check_gko <- function(x, xk){
  # check that xk is a valid knockoff matrix of x
  A <- crossprod(x)
  B <- crossprod(xk)
  # first x^Tx = xk^T xk
  cat("x^Tx = xk^T xk", fill = TRUE)
  cat("Diagonal: ", range(diag(A - B)), fill = TRUE)
  cat("Off-diagonal: ", range(A[upper.tri(A)] - B[upper.tri(B)]), fill = TRUE)

  # then x^Txk = x^tx - S
  C <- crossprod(x, xk)
  cat("x^Txk = x^T x - S", fill = TRUE)
  cat("Off-diagonal: ", range(A[upper.tri(A)] - C[upper.tri(C)]), fill = TRUE)
  cat("Diagonal: ", range(diag(A) - diag(C)), fill = TRUE)
}

debug <- function(){
  load("../conditional/check.RData")
  # we first grab the variables that need to be selected
  xm <- x[, which(group != 0)]
  gm <- group[group != 0]
  # offset: variables that are known to be in the model, and thus will not be selected
  xo <- x[, which(group == 0)]

  # construct group knockoffs
  koequi <- knockoff::create.fixed(X = xm, method = "equi")
  check_gko(koequi$X, koequi$Xk)
  kosdp <- knockoff::create.fixed(X = xm, method = "sdp")
  check_gko(kosdp$X, kosdp$Xk)

  sdp <- function(x) knockoff::create.fixed(x, method = "sdp")
  ko <- knockoff::knockoff.filter(X = xm, y = y, knockoffs = sdp)

  gko <- create_fix_gko(x = xm, y = y, group = gm, method = "equi", sigma = sigma, randomize = FALSE)
  check_gko(gko$x, gko$xk)
}
