
loadings_update <- function(f){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L, D)
    Omega <- (tau * I_n) - (tau^2 * L) %*% Kinv %*% t(L)
    OmegaLD <- Omega %*% L %*% D
    gradL <- -p * (Omega %*% (S %*% OmegaLD) - OmegaLD)   
    
    # line search
    t <- step_size
    line_search <- TRUE
    loss_old <- compute_loss(tau, S, I_n, L, D)
    while(line_search){
      
      # generalized gradient
      G <- (L - prox_fn(L - t * gradL)) / t
      
      # compute losses      
      loss_new <- compute_loss(tau, S, I_n, L - t * G, D)
      
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
    
    rm(Kinv)
    rm(Omega)
    rm(OmegaLD)
    rm(gradL)
    rm(line_search)
    rm(loss_old)
    rm(loss_new)
    rm(G)
    rm(tr)
    rm(l2)
    rm(crt)
    
  }))
}

prior_variance_update <- function(f, q){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L[, -q], D[-q, -q])
    Einv <- tau * I_n - tau^2 * L[, -q] %*% Kinv %*% t(L[, -q])
    
    l <- L[, q]
    u <- Einv %*% l
    
    # names from Tipping and Faul 2003
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

residual_precision_update <- function(f, upper=10.0){
  return(within(f, {
    
    opt <- optim(tau, compute_loss, 
                 method="Brent", 
                 lower=eps, 
                 upper=upper,
                 S=S, I=I_n, L=L, D=D)
    tau <- opt$par
    
    if(opt$convergence != 0){
      stop("residual variance update did not converge")        
    }
    
    rm(opt)
    
  }))
}

mean_update <- function(f){
  return(within(f, {
    
    Kinv <- comp_kinv(tau, L, D)
    Omega <- (tau * I_n) - (tau^2 * L) %*% Kinv %*% t(L)
    Omega_ones <- Omega %*% ones_n
    
    mu <- (t(Y) %*% Omega_ones) / (t(ones_n) %*% Omega_ones)
    S <- comp_samp_cov_wmean(YYt, Y, mu, ones_n)
    
  }))
}
