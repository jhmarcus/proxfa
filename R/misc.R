
comp_samp_cov_wmean <- function(YYt, Y, mu, ones_n){
  
  Y_mu_ones_t <- Y %*% mu %*% t(ones_n)
  S <- YYt - Y_mu_ones_t - t(Y_mu_ones_t) + ones_n %*% t(mu) %*% mu %*% t(ones_n)
  
  return(S/p)
  
}

solve_diag <- function(A, eps=1e-8){
  
  a <- pmax(eps, diag(A))
  Ainv <- diag(1 / a)
  
  return(Ainv)
  
}

tr <- function(A){
  
  return(sum(diag(A)))
  
}

is.sorted <- function(x, ...) {
  
  !is.unsorted(x, ...) | !is.unsorted(rev(x), ...)
  
}


add_convergence_info <- function(f, loss, i){
  
  f$loss <- loss[2:i]
  f$converged <- TRUE
  f$is_decreasing <- is.sorted(f$loss)
  
  return(f)
  
}
