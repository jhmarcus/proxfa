

# used for computing low-rank inverse
comp_kinv <- function(tau, L, D){
  
  Dinv <- solve_diag(D)
  Kinv <- solve(Dinv + tau * t(L) %*% L)
  
  return(Kinv)
  
}


# compute the negative log-likelihood taking avantage of low-rank structure
comp_neg_loglik <- function(tau, S, I, L, D){
  
  Kinv <- comp_kinv(tau, L, D)
  KinvLt <- Kinv %*% t(L)
  SOmega <- (tau * S) - (tau^2) * (S %*% L) %*% KinvLt
  Omega <- (tau * I) - (tau^2) * L %*% KinvLt
  
  tr <- tr(SOmega)
  logdet <- log(det(Omega))
  loss <- (p / 2.0) * (tr - logdet)
  
  return(loss)
  
}