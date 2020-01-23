

# update the loadings with a gradient step either
# with a fixed step-size or line-search
loadings_update <- function(f){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L, D)
    Omega <- (tau * I_n) - (tau^2 * L) %*% Kinv %*% t(L)
    OmegaLD <- Omega %*% L %*% D
    gradL <- -p * (Omega %*% (S %*% OmegaLD) - OmegaLD)   
    
    if(fix_step_size==TRUE){
    
      # gradient-step for fixed step-size  
      L <- prox_fn(L - step_size * gradL)
      
    } else {
      
      # gradient-step using backtracking line-search
      # adapted from: https://www.stat.cmu.edu/~ryantibs/convexopt-S15/lectures/08-prox-grad.pdf
      # TODO: double check 
      t <- step_size
      line_search <- TRUE
      loss_old <- comp_neg_loglik(tau, S, I_n, L, D, p)
      while(line_search){
        
        # generalized gradient
        G <- (L - prox_fn(L - t * gradL)) / t
        
        # compute losses      
        loss_new <- comp_neg_loglik(tau, S, I_n, L - t * G, D, p)
        
        # compute criteria
        tr <- t * tr(t(gradL) %*% G)
        l2 <- .5 * t * sum(G^2)
        crt <- loss_new > (loss_old - tr + l2)
        if(crt){
          
          t <- step_size_shrink * t
          
        } else {
          
          L <- prox_fn(L - t * gradL)
          line_search <- FALSE
          
        }
      }
      
      rm(line_search)
      rm(loss_old)
      rm(loss_new)
      rm(G)
      rm(tr)
      rm(l2)
      rm(crt)
      
    }
    
    rm(Kinv)
    rm(Omega)
    rm(OmegaLD)
    rm(gradL)

  }))
}


# update the qth prior variance holding the other factors fixed 
prior_variance_update <- function(f, q){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L[, -q], D[-q, -q])
    Einv <- tau * I_n - tau^2 * L[, -q] %*% Kinv %*% t(L[, -q])
    
    l <- L[, q]
    u <- Einv %*% l
    
    # naming adopted from Tipping and Faul 2003
    sparsity <- t(l) %*% u
    quality <- tr((S %*% u) %*% t(u))
    if(quality >= sparsity){
      
      D[q, q] <- (quality - sparsity) / (sparsity^2)
      
    } else {
      
      D[q, q] <- eps
      
    }
    
    rm(l)
    rm(Kinv)
    rm(Einv)
    rm(u)
    rm(sparsity)
    rm(quality)
    
  }))
}


# update the residual precision using Brents method
residual_precision_update <- function(f, upper=20.0){
  return(within(f, {
    
    opt <- optim(par=tau, 
                 fn=comp_neg_loglik, 
                 method="Brent", 
                 lower=eps, 
                 upper=upper,
                 S=S, I=I_n, L=L, D=D, p=p)
    
    tau <- opt$par
    
    if(opt$convergence != 0){
      stop("residual variance update did not converge")        
    }
    
    rm(opt)
    
  }))
}


# update the mean for all the features at once and then recompute the 
# sample covariance matrix after removing the mean taking advatage of 
# low rank structure 
mean_update <- function(f){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L, D)
    Omega <- (tau * I_n) - (tau^2 * L) %*% Kinv %*% t(L)
    Omega_ones <- Omega %*% ones_n
    mu <- as.vector(t(Y) %*% Omega_ones / drop(t(ones_n) %*% Omega_ones))
    S <- comp_samp_cov_wmean(YYt, Y, mu, ones_n)
    
  }))
}
