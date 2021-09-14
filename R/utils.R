with_seed <- function(seed, code) {
  code <- substitute(code)
  orig.seed <- .Random.seed
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed)
  eval.parent(code)
}

standardize <- function(x, center = TRUE){
  x.centered <- scale(x, center = center, scale = F)
  x.scaled <- scale(x.centered, center = F, scale = sqrt(colSums(x.centered^2)))
  x.scaled[, ]
}

decomp <- function(x, randomize = FALSE){
  n = nrow(x)
  p = ncol(x)
  stopifnot(n >= 2 * p)
  result = canonical_svd(x)
  Q = qr.Q(qr(cbind(result$u, matrix(0, n, p))))
  u_perp = Q[, (p + 1):(2 * p)]
  if (randomize) {
    Q = qr.Q(qr(rnorm_matrix(p, p)))
    u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}

canonical_svd <- function(x){
  x.svd = tryCatch({
    svd(x)
  }, warning = function(w) {
  }, error = function(e) {
    stop("SVD failed in the creation of fixed-design knockoffs. Try upgrading R to version >= 3.3.0")
  }, finally = {
  })
  for (j in 1:min(dim(x))) {
    i = which.max(abs(x.svd$u[, j]))
    if (x.svd$u[i, j] < 0) {
      x.svd$u[, j] = -x.svd$u[, j]
      x.svd$v[, j] = -x.svd$v[, j]
    }
  }
  return(x.svd)
}
