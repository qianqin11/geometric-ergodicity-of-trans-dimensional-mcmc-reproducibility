###### Mixture model example #####
# Results found in Section 4.3.2.

# log likelihood
# identity is denoted by A in the manuscript. mu is denoted by U in the manuscript.
loglikelihood <- function(y,identity,mu,tau,w) {
  # identity is n dimensional;  indicates the cluster that point i belongs to
  n <- length(y)
  k <- length(mu)
  value <- 0
  nj <- sapply(seq(1,k), function(j){sum(identity==j)}) # k dimensional; indicates the number of points in cluster j
  value <- value + sum(nj*log(w))
  value <- value + sum(dnorm(y, mu[identity], sd=sqrt(tau[identity]), log=T))
  value
}

# log prior
# In the manuscript, mu0 is u_0, a0 is c, b0 is b, alpha0 is gamma.
log.prior <- function(mu,tau,w,k,mu0,tau0,a0,b0,alpha0,pk) {
  value <- log(pk[k]) + lgamma(k+1) + sum(dnorm(mu,mu0,sqrt(tau0*tau), log=T)) +
    k*a0*log(b0) - k*lgamma(a0) - (a0+1)*sum(log(tau)) - (b0)*sum(1/tau) + 
    lgamma(k*alpha0) - k*lgamma(alpha0) + (alpha0-1)*sum(log(w))
  value
}

# \tilde{\pi}_A is the product of logpid.j
logpid.j <- function(nj,sj,ssj,n, mu0, tau0, a0, b0, alpha0) {
  lgamma(nj+alpha0) - (1/2)*log(tau0*nj+1) + 
    lgamma(nj/2+a0) - (nj/2+a0)*log( ssj/2 + mu0^2/2/tau0 -
                                       (sj+mu0/tau0)^2/2/(nj+1/tau0) + b0 )
}



# mixRJMCMC: function for simulating the reversible jump chain for the Gaussian mixture model.
mixRJMCMC <- function(chainlength, y, mu0, tau0, a0, b0, alpha0, pk, 
                      k.ini, id.ini, w.ini, tau.ini, mu.ini, 
                      jumppar=1/3) {
  
  kmax <- length(pk)
  n <- length(y)
  
  idchain <- matrix(rep(0,chainlength*n), nrow=chainlength) # A(t)
  muchain <- matrix(rep(0,chainlength*kmax), nrow=chainlength) # U(t)
  tauchain <- matrix(rep(1,chainlength*kmax), nrow=chainlength) #T(t)
  wchain <- muchain # W(t)
  kchain <- rep(1,chainlength) #K(t)
  
  # initialization:
  idchain[1,] <- id.ini
  muchain[1,1:k.ini] <- mu.ini
  tauchain[1,1:k.ini] <- tau.ini
  wchain[1,1:k.ini] <- w.ini
  kchain[1] <- k.ini
  
  # simulation:
  for (t in 1:(chainlength-1)) {

    if (kchain[t] == 1) {
      qB <- jumppar*min(1, pk[2]/pk[1])
      qD <- 0
    } else if (kchain[t] == kmax) {
      qB <- 0
      qD <- jumppar*min(1, pk[kmax-1]/pk[kmax])
    } else {
      qB <- jumppar*min(1, pk[kchain[t]+1]/pk[kchain[t]])
      qD <- jumppar*min(1, pk[kchain[t]-1]/pk[kchain[t]])
    }
    
    jumprv <- runif(1,0,1)
    
    if (jumprv < qB) {
      # birth move
      k.current <- kchain[t]
      id.current <- idchain[t,]
      w.current <- wchain[t,1:k.current]
      tau.current <- tauchain[t,1:k.current]
      mu.current <- muchain[t,1:k.current]
      
      
      # number of empty clusters:
      k.empty <- k.current - length(unique(id.current))
      
      
      w.birth <- rbeta(1,1,n)
      tau.birth <- 1/rgamma(1,a0,b0)
      mu.birth <- rnorm(1,mu0,sqrt(tau0*tau.birth))
      
      
      k.pro <- k.current+1
      
      arate <- min(1, exp(
        log(k.pro) + lgamma(k.pro*alpha0) - lgamma(alpha0) - lgamma(k.current*alpha0) +
          (n + k.current*alpha0-1)*log(1 - w.birth) + (alpha0-1)*log(w.birth) -
          log(k.empty+1) - dbeta(w.birth,1,n, log=T)
      ))
      
      
      u <- runif(1,0,1)
      
      if(u < arate) {
        id.pro <- id.current
        w.pro <- c((1-w.birth)*w.current, w.birth)
        tau.pro <- c(tau.current,tau.birth)
        mu.pro <- c(mu.current,mu.birth)
        
        l <- 1 # order of mu.pro
        
        
        for (j in 1:k.current) {
          if (mu.birth > mu.current[j])
            l <- j+1
        }
        
        sigma <- 1:k.pro
        sigma[l] <- k.pro
        if (l < k.pro) sigma[(l+1):k.pro] <- l:k.current
        sigma.inv <- 1:k.pro
        sigma.inv[k.pro] <- l
        sigma.inv[l:k.current] <- (l+1):k.pro
        
        id.pro <- sigma.inv[id.pro]
        w.pro <- w.pro[sigma]
        tau.pro <- tau.pro[sigma]
        mu.pro <- mu.pro[sigma]
        
        k.current <- k.pro
        id.current <- id.pro
        w.current <- w.pro
        tau.current <- tau.pro
        mu.current <- mu.pro
      }
      
      
      #####
      
      kchain[t+1] <- k.current
      idchain[t+1,] <- id.current
      wchain[t+1,1:k.current] <- w.current
      tauchain[t+1,1:k.current] <- tau.current
      muchain[t+1,1:k.current] <- mu.current
      
      
    } else if (jumprv >= qB & jumprv < qB+qD) {
      # death move
      
      k.current <- kchain[t]
      id.current <- idchain[t,]
      w.current <- wchain[t,1:k.current]
      tau.current <- tauchain[t,1:k.current]
      mu.current <- muchain[t,1:k.current]
      
      isempty <- rep(T,k.current)
      
      for (i in 1:n) {
        isempty[id.current[i]] <- F
      }
      k.empties <- which(isempty==T)
      k.empty <- length(k.empties)
      
      
      if (k.empty > 0) {
        if (k.empty == 1) {l <- k.empties} else{
          l <- sample(k.empties,1)
        }
        w.death <- w.current[l]
        
        k.pro <- k.current - 1
        
        
        arate <- min(1, exp(
          -log(k.current) - lgamma(k.current*alpha0) + lgamma(alpha0) + lgamma(k.pro*alpha0) -
            (n + k.pro*alpha0-1)*log(1 - w.death) - (alpha0-1)*log(w.death) +
            log(k.empty) + dbeta(w.death,1,n,log=T)
        ))
        
        u <- runif(1,0,1)
        
        if (u < arate) {
          id.pro <- id.current
          if (sum(id.current > l)>0) id.pro[id.current > l] <- id.pro[id.current > l]-1
          w.pro <- w.current[-l]/(1-w.current[l])
          tau.pro <- tau.current[-l]
          mu.pro <- mu.current[-l]
          
          k.current <- k.pro
          id.current <- id.pro
          w.current <- w.pro
          tau.current <- tau.pro
          mu.current <- mu.pro
        }
      }
      
      kchain[t+1] <- k.current
      idchain[t+1,] <- id.current
      wchain[t+1,1:k.current] <- w.current
      tauchain[t+1,1:k.current] <- tau.current
      muchain[t+1,1:k.current] <- mu.current
      
      
    } else {
      # update move
      
      
      sigma <- sample(seq(1,kchain[t]))
      id <- sigma[idchain[t,]]
      id <- idchain[t,]
      k <- kchain[t]
      
      nj <- sapply(seq(1,k), function(j) { sum(id==j) })
      sj <- sapply(seq(1,k), function(j) { sum(y[id==j]) })
      ssj <- sapply(seq(1,k), function(j) { sum(y[id==j]^2) })
      
      
      for (i in 1:n) {  
        # propose new id for a point
        j.pro <- sample(seq(1,kchain[t]), 1)
        
        if (id[i] != j.pro) {
          nj.pro <- nj
          nj.pro[id[i]] <- nj.pro[id[i]]-1
          nj.pro[j.pro] <- nj.pro[j.pro]+1
          sj.pro <- sj
          sj.pro[id[i]] <- sj[id[i]] - y[i]
          sj.pro[j.pro] <- sj[j.pro] + y[i]
          ssj.pro <- ssj
          ssj.pro[id[i]] <- ssj[id[i]] - y[i]^2
          ssj.pro[j.pro] <- ssj[j.pro] + y[i]^2
          
          arate <- min(1, exp(
            logpid.j(nj.pro[j.pro],sj.pro[j.pro],ssj.pro[j.pro],n, mu0, tau0, a0, b0, alpha0) -
              logpid.j(nj[j.pro],sj[j.pro],ssj[j.pro],n, mu0, tau0, a0, b0, alpha0) +
              logpid.j(nj.pro[id[i]],sj.pro[id[i]],ssj.pro[id[i]], n, mu0, tau0, a0, b0, alpha0) -
              logpid.j(nj[id[i]],sj[id[i]],ssj[id[i]], n, mu0, tau0, a0, b0, alpha0)
          ))
          
          u <- runif(1,0,1)
          
          if (u < arate) {
            id[i] <- j.pro
            nj <- nj.pro
            sj <- sj.pro
            ssj <- ssj.pro
          }
        }
      }
      
      
      
      # update w
      w <- rdirichlet(1, nj+alpha0)
      
      # update tau
      tau <- 1/rgamma(kchain[t], a0+nj/2,
                      ssj/2+mu0^2/2/tau0-
                        (sj+mu0/tau0)^2/2/(nj+1/tau0) + b0)
      
      # update mu
      mu <- rnorm(kchain[t], (sj+mu0/tau0)/(nj+1/tau0), sd=sqrt(tau/(nj+1/tau0)))
      
      # re-order
      sigma <- order(mu)
      sigma.inv <- sigma
      sigma.inv[sigma] <- seq(1,kchain[t])
      
      kchain[t+1] <- kchain[t]
      idchain[t+1,] <- sigma.inv[id]
      wchain[t+1,1:kchain[t]] <- w[sigma]
      tauchain[t+1,1:kchain[t]] <- tau[sigma]
      muchain[t+1,1:kchain[t]] <- mu[sigma]
      
      
    }
    
  }
  
  list(k=kchain,id=idchain,w=wchain,tau=tauchain,mu=muchain)
  
}
