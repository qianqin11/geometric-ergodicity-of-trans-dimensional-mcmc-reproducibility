##### functions for simulations in Section 4.2 #######

# If we add a covariate x and its coefficient bscalar to the model, whose X\beta is Xb.pre,
# how much does the log likelihood change?

loglikelihood <- function(x,bscalar, Xb.pre, y) {
  
  fmeans <- as.vector(Xb.pre + cbind(x)%*%cbind(bscalar))
  
  logli <- sum(y*pnorm(fmeans, 0, 1, log.p = TRUE) ) + 
    sum((1-y)*pnorm(fmeans, 0, 1, log.p =TRUE, lower.tail = FALSE) )
  
  logli
}

# The derivative of the log likelihood with respect to bscalar.

dloglikelihood <- function(x, bscalar, Xb.pre, y) {
  fmeans <- as.vector(Xb.pre + cbind(x)%*%cbind(bscalar))
  
  sum(y*exp(dnorm(fmeans, log = TRUE) - pnorm(fmeans, log.p = TRUE))*x) - 
    sum((1-y)*exp(dnorm(fmeans, log = TRUE) - pnorm(fmeans, log.p = TRUE, lower.tail = FALSE))*x)
}

# Second order derivative

ddloglikelihood <- function(x, bscalar, Xb.pre, y) {
  fmeans <- as.vector(Xb.pre + cbind(x)%*%cbind(bscalar))
  fF1 <- exp(dnorm(fmeans, log = TRUE) - pnorm(fmeans,log.p = TRUE))
  fF0 <- exp(dnorm(fmeans, log = TRUE) - pnorm(fmeans,log.p = TRUE, lower.tail=FALSE))
  - sum( y*(fmeans*fF1 + fF1^2)*x^2 ) - sum((1-y)*(-fmeans*fF0+fF0^2)*x^2)
}


# probit.mcmc: generate a reveresible jump Markov chain for the Bayesian probit regression.
# jumppar: 
probit.mcmc <- function(chainlength, k.ini, a.ini, b.ini,
                        X, y, 
                        p, sigma, jumppar=1/3) {
  n <- nrow(X)
  r <- ncol(X)-1
  
  #jump probabilities
  qU <- rep(0,r+1)
  qB <- qU
  qD <- qU
  for (ik in 0:r) {
    qB[ik+1] <- jumppar*min(1, p*(r-ik)/(ik+1))
    qD[ik+1] <- jumppar*min(1, ik/p/(r-ik+1))
    qU[ik+1] <- 1 - qB[ik+1] - qD[ik+1]
  } 
  
  # list of outputs
  models <- matrix(rep(0, (chainlength+1)*r), nrow = chainlength+1)
  intercepts <- rep(0, chainlength+1)
  bs <- matrix(rep(0, (chainlength+1)*r), nrow = chainlength+1)
  jumpproposals <- 0
  jumpacceptances <- 0
  
  # auxillaries
  Xb <- X %*% rbind(intercepts[1], cbind(bs[1,]))
  
  
  # initialization
  models[1,] <- k.ini
  intercepts[1] <- a.ini
  bs[1,] <- as.vector(b.ini)
  
  # MCMC algorithm
  for (t in 2:(chainlength+1)) {
    
    k.current <- models[t-1,]
    Ik.current <- sum(k.current)
    a.current <- intercepts[t-1]
    b.current <- bs[t-1,]
    Xk <- X[,c(1,which(k.current==1)+1)]
    
    move <- runif(1,0,1)
    
    
    if (move < qU[Ik.current+1]) {
      # Update move
      umeans <- as.vector(Xb)
      lowerlimits <- rep(0, n)
      lowerlimits[which(y==0)] <- -Inf
      upperlimits <- rep(Inf, n)
      upperlimits[which(y==0)] <- 0
      u <- rtruncnorm(n, lowerlimits, upperlimits, umeans, sd=1)
      
      zcov <- solve(t(Xk)%*%Xk + 1/sigma^2*diag(Ik.current+1))
      zmean <- zcov %*% t(Xk) %*% u
      z.new <- as.vector(rmvnorm(1, zmean, zcov))
      
      models[t,] <- k.current
      intercepts[t] <- z.new[1]
      bs[t, which(k.current==1)] <- z.new[2:(Ik.current+1)]
      
      Xb <- Xk %*% cbind(z.new)
      
    } else if (move < qU[Ik.current+1] + qB[Ik.current+1]) {
      # birth move
      jumpproposals <- jumpproposals + 1
      if (length(which(k.current == 0))==1) {j <- which(k.current == 0)} else {
        j <- sample(which(k.current == 0), 1)
      }
      
      k.proposed <- k.current
      k.proposed[j] <- 1
      n.logpost <- function(bscalar) {
        - sum(k.proposed) * (log(p) - sqrt(2*pi) - log(sigma) ) + a.current^2/2/sigma^2 + 
          sum(b.current^2)/2/sigma^2 + bscalar^2/2/sigma^2 - 
          loglikelihood(X[,j+1], bscalar, Xb, y)
      }
      n.dlogpost <- function(bscalar) {
        bscalar/sigma^2 - dloglikelihood(X[,j+1], bscalar, Xb, y)
      }
      
      # specify the parameters of the normal proposal of bstar, which is beta.proposed here.
      # For a scalar b, let z(b) be the \beta obtained from changing the jth coordinate of b.current to b.
      # The mean of the normal is the maximizer of b \mapsto \log \pi(k+1,z(b)|y)
      # The variance is the negative second order derivative evaulated at the maximum.
      bscalarmean <- optim(0, n.logpost, n.dlogpost, method = "BFGS")$par
      bscalarvar <- 1/(1/sigma^2-ddloglikelihood(X[,j+1], bscalarmean, Xb, y))
      bscalar.proposed <- rnorm(1, bscalarmean, sd=sqrt(bscalarvar))
      b.proposed <- b.current
      b.proposed[j] <- bscalar.proposed
      
      logpost.current <- sum(k.current) * (log(p) - sqrt(2*pi) - log(sigma) ) - a.current^2/2/sigma^2 - 
        sum(b.current^2)/2/sigma^2 + loglikelihood(X[,j+1], 0, Xb, y)
      logpost.proposed <- -n.logpost(bscalar.proposed)
      acceptprob <- min(1, exp( logpost.proposed +
                                  log(qD[Ik.current+2]) - 
                                  logpost.current -
                                  log(qB[Ik.current+1]) - 
                                  dnorm(bscalar.proposed, bscalarmean, sd=sqrt(bscalarvar), log=TRUE)
      ) / (Ik.current+1) * (r-Ik.current)
      )
      u <- runif(1,0,1)
      if (u < acceptprob) {
        jumpacceptances <- jumpacceptances + 1
        models[t,] <- k.proposed
        intercepts[t] <- a.current
        bs[t,] <- b.proposed
        Xb <- Xb + X[,j+1]%*%cbind(bscalar.proposed)
      } else {
        models[t,] = k.current
        intercepts[t] <- a.current
        bs[t,] <- b.current
      } 
      
    } else {
      # death move
      jumpproposals <- jumpproposals + 1
      if (length(which(k.current == 1)) == 1) {j <- which(k.current == 1)} else {
        j <- sample(which(k.current == 1), 1)
      }
      
      k.proposed <- k.current
      k.proposed[j] <- 0
      bscalar.proposed <- b.current[j]
      b.proposed <- b.current*k.proposed
      Xb.proposed <- Xb - X[,j+1]%*%cbind(bscalar.proposed)
      
      n.logpost <- function(bscalar) {
        - sum(k.current) * (log(p) - sqrt(2*pi) - log(sigma) ) + a.current^2/2/sigma^2 + 
          sum(b.proposed^2)/2/sigma^2 + bscalar^2/2/sigma^2 - 
          loglikelihood(X[,j+1], bscalar, Xb.proposed, y)
      }
      n.dlogpost <- function(bscalar) {
        bscalar/sigma^2 - dloglikelihood(X[,j+1], bscalar, Xb.proposed, y)
      }
      
      
      bscalarmean <- optim(0, n.logpost, n.dlogpost, method="BFGS")$par
      bscalarvar <- 1/(1/sigma^2-ddloglikelihood(X[,j+1], bscalarmean, Xb.proposed, y))
      
      logpost.current <- -n.logpost(bscalar.proposed)
      logpost.proposed <- sum(k.proposed) * (log(p) - sqrt(2*pi) - log(sigma) ) - a.current^2/2/sigma^2 - 
        sum(b.proposed^2)/2/sigma^2 + loglikelihood(X[,j+1], 0, Xb.proposed, y)
      acceptprob <- min(1, exp(
        logpost.proposed +
          log(qB[Ik.current]) +
          dnorm(bscalar.proposed, bscalarmean, sd=sqrt(bscalarvar), log=TRUE) - 
          logpost.current - 
          qD[Ik.current+1]
      ) / (r-Ik.current+1) * Ik.current
      )
      
      u <- runif(1,0,1)
      if (u < acceptprob) {
        jumpacceptances <- jumpacceptances + 1
        models[t,] <- k.proposed
        intercepts[t] <- a.current
        bs[t,] <- b.proposed
        Xb <- Xb - X[,j+1] %*% cbind(bscalar.proposed)
      } else {
        models[t,] = k.current
        intercepts[t] <- a.current
        bs[t,] <- b.current
      }    
    }
    
  }
  
  
  list(k=models,a=intercepts,b=bs,
       jumpproposals=jumpproposals,jumpacceptances=jumpacceptances)
}




