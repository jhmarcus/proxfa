
box_prox_fn <- function(A, lower=1e-8, upper=1-1e-8){
  
  A[A<=lower] <- lower
  A[A>=upper] <- upper
  
  return(A)
  
}


non_negative_prox_fn <- function(A, lower=1e-8){
  
  A[A<=lower] <- lower

  return(A)
}


row_simplex_prox_fn <- function(A){
  
  n <- nrow(A)
  for(j in 1:n){
    A[j, ] <- alstructure:::projsplx(A[j, ])   
  }
 
  return(A)
  
}


col_simplex_prox_fn <- function(A){
  
  n <- ncol(A)
  for(j in 1:n){
    A[, j] <- alstructure:::projsplx(A[, j])   
  }
  
  return(A)
  
}


identity_prox_fn <- function(A){
  
  return(A)
  
}