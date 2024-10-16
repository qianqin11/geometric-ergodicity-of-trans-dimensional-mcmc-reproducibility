##### Robust autoregression example ######
# Results found in Section d of Supplement II

# log posterior density
lpost.density <- function(X.full, beta.full, y, u, tau, a.pr, sigma.pr) {
  n <- length(y)
  q <- length(beta.full)
  logli <- -n/2* log(tau) - (3/2)*sum(log(u)) - (1/8)*sum(1/u) - 
    (1/2/tau)*sum(u*(y-as.vector(X.full%*%cbind(beta.full)))^2)
  logpr <- -(a.pr+1)/2*log(tau) - q/2*log(2*pi) - q*log(sigma.pr) - q/2*log(tau)
  - (1/2/sigma.pr^2/tau)*sum(beta.full^2)
  logli + logpr
}

# Function for generating the chain
# ini means initial values.
# prior.k is the prior density for K
# prior density for \tau is tau^(-(prior.a+1)/2)
# prior density for \mu given \tau is N(0,\tau*prior.scale)
# 1-2*jumppar: lowest probability of conducting an update move.
robust.mcmc <- function(chainlength, k.ini, a.ini, b.ini, tau.ini,
                        X, y, y.start, prior.k, prior.a = 1, prior.scale,
                        jumppar = 1/3, silent = FALSE
) {
  n <- nrow(X)
  p <- ncol(X)
  kmax <- length(prior.k)-1
  start.length <- length(y.start)
  y.start <- y.start[(start.length - kmax+1):start.length]
  y.full <- c(y.start, y)
  ymat <- matrix(y, ncol = 1)
  
  
  W <- matrix(rep(0, n*kmax), nrow = n)
  for (i in 1:n) {
    for (k in 1:kmax) {
      W[i, k] <- y.full[i+kmax-k]
    }
  }
  
  X.full <- cbind(X, W)
  
  # jump probabilities
  qU <- rep(0, kmax + 1)
  qB <- rep(0, kmax + 1)
  qD <- rep(0, kmax + 1)
  for (k in 0:kmax) {
    if (k < kmax) {
      qB[k+1] <- jumppar*min(1, prior.k[k+2]/prior.k[k+1])
      qD[k+2] <- jumppar*min(1, prior.k[k+1]/prior.k[k+2])
      qU[k+1] <- 1 - qB[k+1] - qD[k+1]
    } else {
      qB[k+1] <- 0
      qU[k+1] <- 1 - qB[k+1] - qD[k+1]
    }
  }
  
  # list of outputs
  models <- rep(0,chainlength+1)
  avector <- matrix(rep(0, kmax*(chainlength+1)), ncol = kmax)
  bvector <- matrix(rep(0, p*(chainlength+1)), ncol = p)
  uvector <- matrix(rep(0, n*(chainlength+1)), ncol = n)
  tauvector <- rep(0, chainlength)
  
  jumpproposals <- 0
  acceptances <- 0
  
  # initialization
  models[1] <- k.ini
  if (k.ini >= 1) avector[1,1:k.ini] <- a.ini
  bvector[1,] <- b.ini
  tauvector[1] <- tau.ini
  
  
  
  # MCMC algorithm
  for (t in 2:(chainlength+1)) {
    
    if (silent==FALSE) {
      if ((t-1)%%1000 == 0) print(t-1)
    }
    
    k.current <- models[t-1]
    Xfull.current <- X.full[, 1:(p+k.current), drop=FALSE]
    b.current <- bvector[t-1,]
    if (k.current >= 1) {
      a.current <- avector[t-1,1:k.current]
      bfull.current <- c(b.current,a.current)
    } else {
      bfull.current <- b.current
    }
    
    tau.current <- tauvector[t-1]
    u.current <- uvector[t-1,]
    
    move <- runif(1,0,1)
    if (move < qU[k.current+1]) {
      # update move
      
      
      u.new <- as.vector(rinvgauss(n, mean = sqrt(tau.current)/2/abs(y - as.vector(Xfull.current%*%cbind(bfull.current))),
                                   dispersion = 4))
      
      
      Q <- diag(u.new, ncol = n)
      XQy <- t(Xfull.current) %*% Q %*% ymat
      XQX.inv <- solve(t(Xfull.current)%*%Q%*%Xfull.current + diag(1/prior.scale^2, nrow = p+k.current) )
      
      tau.new <- 1/rgamma(1, (n+prior.a-1)/2, 
                          rate = (1/2)*as.numeric( t(ymat)%*%Q%*%ymat - 
                                                     t(XQy)%*%XQX.inv%*%XQy  )  )
      
      if (is.nan(tau.new)) tau.new <- 1e100
      
      bfull.new <- as.vector(rmvnorm(1, XQX.inv%*%XQy, tau.new*XQX.inv))
      
      models[t] <- k.current
      if(k.current >= 1) {avector[t, 1:k.current] <- bfull.new[(p+1):(p+k.current)]}
      bvector[t,] <- bfull.new[1:p]
      uvector[t,] <- u.new
      tauvector[t] <- tau.new
      
      
    } else if (move < qU[k.current+1] + qB[k.current+1]) {
      # birth move
      
      jumpproposals <- jumpproposals + 1
      
      k.proposed <- k.current + 1
      wkplusone <- W[,k.proposed]
      wkplusone2 <- sum(u.current*wkplusone^2)
      
      amean <- (1/(wkplusone2+1/prior.scale^2))*as.numeric( rbind(wkplusone*u.current) %*% (ymat - Xfull.current %*% cbind(bfull.current) )  )
      asd <- sqrt(1/(wkplusone2+1/prior.scale^2))
      
      akplusone <- rnorm(1, amean, sd = asd)
      
      Xfull.proposed <- X.full[, 1:(p+k.proposed)]
      bfull.proposed <- c(bfull.current, akplusone)
      
      lpost.current <- lpost.density(Xfull.current, bfull.current, y, u.current, tau.current, prior.a, prior.scale)
      lpost.proposed <- lpost.density(Xfull.proposed, bfull.proposed, y, u.current, tau.current, prior.a, prior.scale)
      
      acceptprob <- min(1, exp(lpost.proposed-lpost.current-
                                 dnorm(akplusone, amean, asd, log = TRUE)))
      
      if(!is.numeric(acceptprob) | is.nan(acceptprob)) acceptprob <- 0
      
      rand <- runif(1,0,1)
      if (rand < acceptprob) {
        models[t] <- k.proposed
        avector[t,] <- avector[t-1,]
        avector[t, k.proposed] <- akplusone
        bvector[t,] <- b.current
        uvector[t,] <- u.current
        tauvector[t] <- tau.current
        acceptances <- acceptances + 1
      } else {
        models[t] <- k.current
        avector[t, ] <- avector[t-1, ]
        bvector[t,] <- b.current
        uvector[t,] <- u.current
        tauvector[t] <- tau.current
      }
      
      
    } else {
      # death move
      
      jumpproposals <- jumpproposals + 1
      
      k.proposed <- k.current - 1
      wk <- W[, k.current]
      wk2 <- sum(wk^2*u.current)
      
      Xfull.proposed <- X.full[, 1:(p+k.proposed)]
      bfull.proposed <- bfull.current[1:(p+k.proposed)]
      
      amean <- (1/(wk2+1/prior.scale^2))*as.numeric( rbind(wk*u.current) %*% (ymat - Xfull.proposed %*% cbind(bfull.proposed) )  )
      asd <- sqrt(1/(wk2+1/prior.scale^2))
      
      #ak <- a.current[k.proposed]
      ak <- bfull.current[p+k.current]
      
      lpost.current <- lpost.density(Xfull.current, bfull.current, y, u.current, tau.current, prior.a, prior.scale)
      lpost.proposed <- lpost.density(Xfull.proposed, bfull.proposed, y, u.current, tau.current, prior.a, prior.scale)
      
      acceptprob <- min(1, exp(lpost.proposed-lpost.current + dnorm(ak, amean, asd, log = TRUE)))
      if(!is.numeric(acceptprob) | is.nan(acceptprob)) acceptprob <- 0
      
      rand <- runif(1,0,1)
      if (rand < acceptprob) {
        models[t] <- k.proposed
        if (k.proposed >= 1) {
          avector[t, 1:k.proposed] <- a.current[1:k.proposed]
        } 
        bvector[t,] <- b.current
        uvector[t,] <- u.current
        tauvector[t] <- tau.current
        acceptances <- acceptances + 1
      } else {
        models[t] <- k.current
        avector[t, 1:k.current] <- a.current
        bvector[t,] <- b.current
        uvector[t,] <- u.current
        tauvector[t] <- tau.current
      }
      
    }
    
  } 
  
  list(k=models,b=bvector,a=avector,t=tauvector,
       jumpproposals=jumpproposals,jumpacceptances=acceptances)
  
}
