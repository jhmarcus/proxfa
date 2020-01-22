

# fast computation of sample covariance matrix with feature
# specific mean without dense matrix-matrix multiplication
comp_samp_cov_wmean <- function(YYt, Y, mu, ones_n){
  
  Y_mu_ones_t <- Y %*% mu %*% t(ones_n)
  S <- YYt - Y_mu_ones_t - t(Y_mu_ones_t) + ones_n %*% t(mu) %*% mu %*% t(ones_n)
  
  return(S/p)
  
}


# invert a diagonal matrix
solve_diag <- function(A, eps=1e-8){
  
  a <- diag(A)
  a[a<=eps] <- eps
  Ainv <- diag(1 / a)
  
  return(Ainv)
  
}


# compute the trace of a matrix A
tr <- function(A){
  
  return(sum(diag(A)))
  
}


# check if a vector is sorted from:
# https://stackoverflow.com/questions/23547929/determine-if-a-vector-is-ordered
is.sorted <- function(x, ...) {
  
  !is.unsorted(x, ...) | !is.unsorted(rev(x), ...)
  
}


# adds some additional information like the loss to 
# the return list
add_convergence_info <- function(f, loss, i){
  
  f$loss <- loss[2:i]
  f$converged <- TRUE
  f$is_decreasing <- is.sorted(f$loss)
  f$i <- i
  
  return(f)
  
}