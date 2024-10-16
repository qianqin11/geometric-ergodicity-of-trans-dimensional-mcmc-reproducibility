### This files contains functions for uncertainty assessment (confidence interval construction).

# batch means estimator
batchmeans <- function(chain, bn=floor(nrow(as.matrix(chain))^(0.6))) {
  # chainmat stores the chain. Each row of chainmat corresponds to one iteration.
  chainmat <- as.matrix(chain)
  p <- ncol(chainmat)
  n <- nrow(chainmat)-1
  an <- n%/%bn
  chain.truncated <- chainmat[(n-an*bn+2):(n+1),,drop=FALSE]
  globalmean <- apply(X=chain.truncated, MARGIN = 2, FUN=mean)
  asy.var <- rep(0,p)
  for (i in 1:an) {
    localmean <- apply(X=chain.truncated[((i-1)*bn+1):(i*bn),,drop=FALSE], MARGIN = 2, FUN = mean)
    asy.var <- asy.var + (localmean-globalmean)^2
  }
  asy.var <- asy.var*bn/(an-1)
  c(asy.var,bn)
}


# Construct Wald intervals
# inflate means adding (\log m) \sqrt(1/b_m^2 + b_m/m) to the batch-means estimator.
# adjusted means using Bonferroni adjustment.
get.intervals <- function(alpha, samplemean, vbatch, n, inflate=T, silent = TRUE, adjusted=F) {
  m <- length(samplemean)
  
  if (inflate==T) {
    v <- vbatch[1:m] + log(n)*sqrt(1/vbatch[m+1]^2+vbatch[m+1]/n)
  } else {
    v <- vbatch[1:m]
  }
  
  if (adjusted == T) {
    z <- qnorm(1-alpha/2/m)
  } else {
    z <- qnorm(1-alpha/2)
  }
  
  
  lowers <- samplemean-z*sqrt(v)/sqrt(n)
  uppers <- samplemean+z*sqrt(v)/sqrt(n)
  matrix(c(lowers, samplemean, uppers), ncol=3)
}