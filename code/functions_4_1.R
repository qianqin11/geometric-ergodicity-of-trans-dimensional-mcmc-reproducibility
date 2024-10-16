### This file contains functions that are used for the simulations in Section 4.1.

# norm returns the norm of a Markov transition matrix (Mtm) P whose stationary distribution is pivector.
norm <- function(P,pivector=rep(1,nrow(P))) {
  # P is the transition matrix, pivector is the vector of the stationary distribution.
  if (length(P) == 1) return(0)
  D.sqrt <- diag(sqrt(pivector))
  D.sqrt.inv <- diag(1/sqrt(pivector))
  sqrt(eigen(D.sqrt.inv %*% t(P) %*% D.sqrt %*% D.sqrt %*% P %*% D.sqrt.inv )$value[2])
}

# gap returns 1 minus the second largest eigenvalue of an Mtm.
gap <- function(P) {
  if (length(P) == 1) return(1)
  1-eigen(P)$value[2]
}


# generate.mtm: function for generating the Mtm P in Section 4.1.
# Gl: "fast" or "slow", Lo: "fast", "slow" or "varied".
# Zksize is denoted by n in Section 4.1.

generate.mtm <- function(kmax, Zksize, Gl, Lo) {
  # Initialization:
  P <- matrix(rep(0,(kmax*Zksize)^2), nrow=kmax*Zksize)
  
  for (k in 1:kmax) {
    for (z in 1:Zksize) {
      # Label (k,z) using a single index:
      index <- (k-1)*Zksize + z
      
      # Transition probabilities for the global moves:
      if (Gl == "fast") {
        for (k1 in 1:kmax) {
          P[index, (k1-1)*Zksize+z] <- P[index, (k1-1)*Zksize+z]+1/2/kmax
        }
      } else if (Gl == "slow") {
        P[index,index] <- P[index,index]+1/4
        if (k==kmax) {
          P[index, z] <- P[index,z] + 1/8
        } else {
          P[index, index+Zksize] <- P[index,index+Zksize] + 1/8
        }
        if (k==1) {
          P[index, (kmax-1)*Zksize+z] <- P[index, (kmax-1)*Zksize+z] + 1/8
        } else {
          P[index, index-Zksize] <- P[index, index-Zksize] + 1/8
        }
      }
      
      # Transition probabilities for the local moves:
      if (Lo == "fast") {
        for (z1 in 1:Zksize) {
          P[index, (k-1)*Zksize+z1] <- P[index, (k-1)*Zksize+z1] + 1/2/Zksize
        }
      } else if (Lo == "slow") {
        P[index,index] <- P[index,index] + 1/4
        if (z==Zksize) {
          P[index, (k-1)*Zksize+1] <- P[index, (k-1)*Zksize+1] + 1/8
        } else {
          P[index, index+1] <- P[index, index+1] + 1/8
        }
        if (z==1) {
          P[index, k*Zksize] <- P[index, k*Zksize] + 1/8
        } else {
          P[index, index-1] <- P[index, index-1] + 1/8
        }
      } else if (Lo == "varied") {
        for (z1 in 1:Zksize) {
          P[index, (k-1)*Zksize+z1] <- P[index, (k-1)*Zksize+z1] + (1/k)*1/2/Zksize
        }
        P[index,index] <- P[index,index] + 1/2*(1-1/k)
      }
      
    }
  }
  P
}

# get.Pbar: calculate \bar{P} from P.
get.Pbar <- function(P,kmax,Zksize) {
  Pbar <- matrix(rep(0,kmax^2), nrow=kmax)
  for (k1 in 1:kmax) {
    for (k2 in 1:kmax) {
      Pbar[k1,k2] <- (1/Zksize)*sum(P[ ((k1-1)*Zksize+1):(k1*Zksize),
                                       ((k2-1)*Zksize+1):(k2*Zksize)])
    }
  }
  Pbar
}

# get.P2bar: calculate \overline{P^*P} from P.
get.P2bar <- function(P,kmax,Zksize) {
  D <- diag(rep(1,kmax*Zksize), nrow=kmax*Zksize)
  D.inv <- D
  get.Pbar(D.inv%*%t(P)%*%D%*%P,kmax,Zksize)
}

# get.Pk: calculate P_k.
get.Pk <- function(k,Zksize,Lo) {
  Pk <- matrix(rep(0,Zksize^2), nrow=Zksize)
  
  for (z1 in 1:Zksize) {
    if (Lo == "fast") {
      Pk[z1,1:Zksize] <- 1/Zksize
    } else if (Lo == "slow") {
      Pk[z1,z1] <- Pk[z1,z1] + 1/2
      if (z1==Zksize) {
        Pk[z1,1] <- Pk[z1,1] + 1/4
      } else {
        Pk[z1,z1+1] <- Pk[z1,z1+1] + 1/4
      }
      if (z1==1) {
        Pk[z1,Zksize] <- Pk[z1,Zksize] + 1/4
      } else {
        Pk[z1,z1-1] <- Pk[z1,z1-1] + 1/4
      }
    } else if (Lo == "varied") {
      Pk[z1,z1] <- 1-1/k
      Pk[z1,1:Zksize] <- Pk[z1,1:Zksize]+1/k/Zksize
    }
    
  }
  
  Pk
}

# obtaining bounds on 1-\|P\|.
# reversible: exploit reversibility and positive semi-definite-ness?

getbound <- function(kmax,Zksize,Gl,Lo,reversible) {
  P <- generate.mtm(kmax,Zksize,Gl,Lo)
  m <- kmax*Zksize
  Pivector <- rep(1/m,m)
  if(reversible) {
    # reversible bound
    P.norm <- 1-gap(P)
    Pbar.gap <- gap(get.Pbar(P,kmax,Zksize))
    Pk.gaps <- rep(0,kmax)
    for (k in 1:kmax) {
      Pk.gaps[k] <- gap(get.Pk(k,Zksize,Lo))
    }
    Pk.gap <- min(Pk.gaps)
    P.bound <- sqrt(1-1/2*Pk.gap*Pbar.gap)
    return(list(truth=1-P.norm, bound=1-P.bound,Pk=Pk.gap,Pbar=Pbar.gap))
  } else {
    # non-reversible bound
    P.norm <- norm(P, rep(1,m))
    P2bar.gap <- gap(get.P2bar(P,kmax,Zksize))
    Pk.norms <- rep(1,kmax)
    for (k in 1:kmax) {
      Pk.norms[k] <- norm(get.Pk(k,Zksize,Lo), rep(1,Zksize))
    }
    Pk.norm <- max(Pk.norms)
    P.bound <- (1-1/4*(1-Pk.norm^2)*P2bar.gap)^(1/4)
    return(list(truth=1-P.norm, bound=1-P.bound,Pk=1-Pk.norm,Pbar=P2bar.gap))
  }
  # truth is 1-\|P\|, bound is a lower bound on it.
  # Pk and Pbar give 1-\|P_k\|_{\Phi_k} and the spectral gap of \overline{P^*P}.
} 