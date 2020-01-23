
library(MASS)

test_that(paste("the low-rank approach to computing the inverse covariance",
                "is close to the brute force approach"),  {
  n <- 40            
  Q = 5
  p <- 500
  L <- matrix(rnorm(n * Q), nrow=n, ncol=Q)
  D <- diag(rep(1, Q))
  I_n <- diag(rep(1, n))
  tau <- 1
  
  Sigma <- L %*% D %*% t(L) + (1 / tau) * I_n
  Omega_brute <- solve(Sigma)
  
  Kinv <- comp_kinv(tau, L, D)
  Omega_low_rank <- (tau * I_n) - (tau^2 * L) %*% Kinv %*% t(L)
  
  expect_equal(Omega_brute, Omega_low_rank, tolerance = 1e-8)
})


test_that(paste("the low-rank approach to computing the negative-log-likelihood",
                "is close to the brute force approach"),  {
  n <- 40            
  Q = 5
  p <- 500
  L <- matrix(rnorm(n * Q), nrow=n, ncol=Q)
  D <- diag(rep(1, Q))
  I <- diag(rep(1, n))
  tau <- 1
  
  Sigma <- L %*% D %*% t(L) + (1 / tau) *  I
  Omega <- solve(Sigma)
  Y <- t(MASS::mvrnorm(n=p, mu=rep(0, n), Sigma=Sigma))
  S <- Y %*% t(Y) / p 
  
  brute_force_nll <- (p / 2.0) * (tr(S %*% Omega) - log(det(Omega)))
  low_rank_nll <- comp_neg_loglik(tau, S, I, L, D, p)
  
  expect_equal(low_rank_nll, brute_force_nll, tolerance = 1e-8)
})
