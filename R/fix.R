create_equivariant <- function(x, gr, randomize){
  # gr is the sorted grouping information in the form of (1, 1, 1, 2, 2, 3, 3, 3, ...)
  x.svd = decomp(x, randomize)
  # sanity check of svd
  # range(x - x.svd$u %*% diag(x.svd$d) %*% t(x.svd$v))

  if (any(x.svd$d <= 1e-05 * max(x.svd$d)))
    stop(paste("Data matrix is rank deficient.", "Equicorrelated knockoffs will have no power."))

  p <- ncol(x)
  # S and D are the matrices corresponding to the same notation as in Dai&Barbar(2016)
  S <- matrix(0, p, p)
  D <- matrix(0, p, p)

  gg <- unique(gr)
  for(j in gg){
    idx <- which(gr == j)
    if(length(idx) > 1){
      # if j-th group has more than 1 variable
      xsub.svd <- canonical_svd(x = x[, idx])
      S[idx, idx] <- xsub.svd$v %*% diag(xsub.svd$d^2) %*% t(xsub.svd$v)
      D[idx, idx] <- xsub.svd$v %*% diag(1 / xsub.svd$d) %*% t(xsub.svd$v)
    }
    else{
      # if j-th group has only 1 variable
      S[idx, idx] <- crossprod(x[, idx])
      D[idx, idx] <- 1 / sqrt(crossprod(x[, idx]))
    }
  }

  # compute the value of gamma
  DsigD <- eigen(crossprod(x %*% D), only.values = TRUE)
  gam <- min(1, 2 * min(DsigD$values))
  # and thus determine the matrix S
  S <- gam * S

  # construct C
  CtC <- 2 * S - S %*% (x.svd$v %*% diag(1 / (x.svd$d^2)) %*% t(x.svd$v)) %*% S
  eig <- eigen(CtC, only.values = TRUE)
  # to avoid numerical issues
  if(any(eig$values < 0)){
    CtC <- CtC + diag(abs(min(eig$values)) + 1e-4, nrow = p, ncol = p)
  }
  else if(min(eig$values) < 1e-8){
    CtC <- CtC + diag(1e-4, nrow = p, ncol = p)
  }
  C <- chol(CtC)

  # now return knockoffs
  xk <- x - x.svd$u %*% diag(1 / x.svd$d) %*% t(x.svd$v) %*% S + x.svd$u_perp %*% C

  # browser()
  # check_gko(x, xk)
  return(xk)
}
