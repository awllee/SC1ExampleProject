# input: matrix
# output: vector of singular values
get.singular.values <- function(mtx) {
  decomposition <- svd(mtx, 0, 0)
  return(decomposition$d)
}

# input: singular values
# output: list of rank and frobenius norm for each truncated svd
get.frobenius.norms <- function(sigmas) {
  k <- 0:length(sigmas)
  ss <- sum(sigmas^2)
  norm <- c(sqrt(ss), sqrt(ss - cumsum(sigmas^2)))
  return(list(k=k, norm=norm))
}

# input: singular values, no. rows, no. columns, target average error
# output: smallest rank with average error <= target
.min.k.average.error <- function(sigmas, m, n, value) {
  errors <- get.frobenius.norms(sigmas)
  vs <- errors$norm / sqrt(m*n)
  k.index <- which.max(vs <= value)
  return(errors$k[k.index])
}

# input: matrix and target average error
# output: smallest rank with average error <= target
.get.rank.avg.error <- function(mtx, avg.error) {
  sigmas <- get.singular.values(mtx)
  dims <- dim(mtx)
  return(.min.k.average.error(sigmas, dims[1], dims[2], avg.error))
}

# input: matrix and rank between 1 and min(dims(mtx))
# output: compressed matrix with decomposition, rank, error, avg.error
#   and (compression) ratio.
.compress.matrix.rank <- function(mtx, rank) {
  stopifnot(rank > 0 && rank <= min(dim(mtx)))
  decomposition <- svd(mtx, rank, rank)
  sigmas <- decomposition$d
  error <- sqrt(sum(sigmas[(rank+1):length(sigmas)]^2))
  avg.error <- error/sqrt(prod(dim(mtx)))
  dims <- dim(mtx)
  ratio <- prod(dims)/(rank*(sum(dims)+1))
  # remove unnecessary singular values
  decomposition$d <- decomposition$d[1:rank]
  out <- list(decomposition=decomposition, rank=rank, error=error, avg.error=avg.error,
              ratio=ratio)
  class(out) = "compressed.matrix"
  return(out)
}

# input: object
# output: TRUE <=> object is a compressed matrix.
is.compressed.matrix <- function(x) {
  return(class(x) == "compressed.matrix")
}

# input: matrix and target (compression) ratio
# output: largest rank with ratio >= target
.get.rank.ratio <- function(mtx, ratio) {
  dims <- dim(mtx)
  m <- dims[1]
  n <- dims[2]
  p <- min(m, n)
  ks <- 1:p
  ratios <- m*n/(ks*(m+n+1))
  rank.plus.1 <- min(which(ratios < ratio))
  if (rank.plus.1 > 1) {
    return(rank.plus.1 - 1)
  } else {
    warning("ratio cannot be achieved, using rank 1")
    return(rank.plus.1)
  }
}

# input: matrix, and one of the three target arguments
# output: compressed matrix satisfying the given target
compress.matrix <- function(mtx, rank=NA, ratio=NA, avg.error=NA) {
  if (sum(!is.na(c(rank,ratio,avg.error))) != 1) {
    stop("exactly one of rank, ratio and avg.error should be specified")
  }
  if (!is.na(ratio)) {
    rank <- .get.rank.ratio(mtx, ratio)
  }
  if (!is.na(avg.error)) {
    rank <- .get.rank.avg.error(mtx, avg.error)
  }
  if (!is.na(rank)) {
    return(.compress.matrix.rank(mtx, rank))
  }
}

# input: compressed matrix
# output: decompressed matrix
decompress.matrix <- function(cmtx) {
  stopifnot(is.compressed.matrix(cmtx))
  U <- cmtx$decomposition$u
  V <- cmtx$decomposition$v
  sigmas <- cmtx$decomposition$d
  k <- cmtx$rank
  return(U %*% diag(sigmas[1:k],k,k) %*% t(V))
}

# input: compressed matrix
# effect: summary is printed
summary.compressed.matrix <- function(cmtx) {
  s <- paste0("Rank = ", cmtx$rank, "\n")
  s <- paste0(s, "Compression ratio = ", cmtx$ratio, "\n")
  s <- paste0(s, "Using the Frobenius norm:", "\n")
  s <- paste0(s, "  Total error = ", cmtx$error, "\n")
  s <- paste0(s, "  Average error = ", cmtx$avg.error, "\n")
  cat(s)
}
