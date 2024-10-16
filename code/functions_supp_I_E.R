### This file contains functions that are used for the simulations in Section E of Supplement I.

# rjfinite: function for generating the toy chain associated with P when Gl="slow", Lo="slow".
rjfinite <- function(kmax, Zksize, chainlength, k.ini, z.ini) {
  ks <- rep(1,chainlength)
  zs <- rep(1,chainlength)
  
  # initial values:
  ks[1] <- k.ini
  zs[1] <- z.ini
  
  for (t in 2:chainlength) {
    r <- runif(1,0,1)
    if (r < 1/2) {
      # update move
      ks[t] <- ks[t-1]
      r1 <- runif(1,0,1)
      if (r1 < 1/2) {
        # stay
        zs[t] <- zs[t-1]
      } else if (r1 < 3/4) {
        # move right
        if (zs[t-1]==Zksize) {
          zs[t] <- 1
        } else {
          zs[t] <- zs[t-1]+1
        }
      } else {
        # move left
        if (zs[t-1]==1) {
          zs[t] <- Zksize
        } else {
          zs[t] <- zs[t-1]-1
        }
      }
    } else if (r>=1/2) {
      # jump move
      zs[t] <- zs[t-1]
      r2 <- runif(1,0,1)
      if (r2 < 1/2) {
        # stay
        ks[t] <- ks[t-1]
      } else if (r2 < 3/4) {
        # move up
        if (ks[t-1]==kmax) {
          ks[t] <- 1
        } else {
          ks[t] <- ks[t-1] + 1
        } 
      } else {
        # move down
        if (ks[t-1]==1) {
          ks[t] <- kmax
        } else {
          ks[t] <- ks[t-1] - 1
        }
      }
    }
  }
  
  list(k=ks,z=zs)
}

# Function for calculating the empirical coverage rate.
get.coverage.rate <- function(len,trialrep,alpha,kmax,Zksize) {
  
  truth <- rep(1,kmax)/kmax
  
  
  alphalength <- length(alpha)
  
  output <- matrix(rep(0,2*alphalength*5), ncol=5)
  
  colnames(output) <- c("truth", "nominal", "inflate", "empirial", "length")
  
  for (inflate in 0:1) {
    for (j in 1:alphalength) {
      index <- (j-1)*2+inflate+1
      output[index,1] <- truth[1]
      output[index,2] <- 1-alpha[j]
      output[index,3] <- inflate
    }
  }
  
  for (i in 1:trialrep) {
    chain.finite <- rjfinite(kmax,Zksize,len,1,1)
    kind <- matrix(rep(0,kmax*len), ncol=kmax)
    for (t in 1:len) {
      for (j in 1:kmax) {
        if (chain.finite$k[t] == j) {
          kind[t,j] <- 1
        }
      }
    }
    kind.mean <- apply(kind,2,mean)
    kind.var <- batchmeans(kind)
    
    for (inflate in 0:1) {
      for (j in 1:alphalength) {
        intervals.finite <- get.intervals(alpha[j],kind.mean,kind.var,len,inflate=inflate, adjusted=T)
        index <- (j-1)*2+inflate+1
        output[index,5] <- output[index,5] + mean((intervals.finite[,3]-intervals.finite[,1])/2)
        if(prod(truth<intervals.finite[,3] & truth>intervals.finite[,1]) == 1) {
          output[index,4] <- output[index,4] + 1
        }
        
      }
    }
    
    print(i)
  }
  
  output[,4:5] <- output[,4:5]/trialrep
  output
}