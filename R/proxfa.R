


#' @title Proximal gradient descent for factor analysis
#'
#' @description IN PROGRESS
#'
#' @param Y is the n x p data matrix with no missing entries
#' 
#' @param L is the intial n x Q loadings matrix 
#' 
#' @param D is the intial Q x Q diagonal matrix of prior variances
#' 
#' @param mu is the intial p x 1 vector of feature specific means
#' 
#' @param tau is the intial the residual precision
#'
#' @param step_size is the intial step_size for proximal 
#'   gradient descent with line search
#' 
#' @param step_size_shrink is a value between 0 and 1 to 
#'   shrink the step-size in each line search inner iteration
#' 
#' @param max_iter is the maximum number of iterations to run
#'   the algorithim
#'
#' @param max_inner_iter is the maximum number of inner iterations for
#'   the inner loop of proximal gradient descent on the loadings
#'   
#' @param tol is the convergence tolerance for the change in the loss
#' 
#' @param inner_tol is the convergence tolerance for the change in loss for 
#'   the inner loop of proximal gradient descent on the loadings
#'   
#' @param n_print is the interval of iterations before printing an update
#'
#' @param prox_fn is a function that implements the proximal operation on L
#'   by default is is the identity
#' 
#' @param eps is a small number for dealing with numerical precision 
#' 
#' @param fix_D is a boolean to update the prior variances
#' 
#' @param fix_mu is a boolean to update the means
#' 
#' @param fix_tau is a boolean to update the residual precision
#' 
#' @return returns a list with fitted values
#' 
#' @export
proxfa <- function(Y, L, D, mu, tau, 
                   step_size, step_size_shrink=.5,
                   max_iter=500, max_inner_iter=20, 
                   tol=1e-3, inner_tol=1e-3, n_print=10, 
                   prox_fn=function(L){L}, eps=1e-8, fix_D=FALSE, 
                   fix_mu=FALSE, fix_tau=FALSE){
  
  n <- nrow(Y) # number of samples
  p <- ncol(Y) # number of features
  n_factors <- ncol(L) # number of factors
  ones_n <- rep(1, n)
  I_n <- diag(ones_n)
  
  # compute sample covarinace including the feature means
  YYt <- Y %*% t(Y)
  S <- comp_samp_cov_wmean(YYt, Y, mu, ones_n)
  
  # intial fit object
  f <- list(Y=Y, YYt=YYt, n=n, p=p, 
            step_size=step_size, step_size_shrink=step_size_shrink,
            S=S, L=L, D=D, mu=mu, tau=tau, ones_n=ones_n, I_n=I_n, 
            prox_fn=prox_fn,  eps=eps)
            
  loss = rep(NA, max_iter+1)
  loss[1] <- compute_loss(f$tau, f$S, f$I_n, f$L, f$D)
  for(i in 2:(max_iter+1)){
    
    ########## loadings update ########## 
    last_inner_loss <- loss[1]
    
    # gradient descent
    for(m in 1:max_inner_iter){
      f <- loadings_update(f)  
      inner_loss <- compute_loss(f$tau, f$S, f$I_n, f$L, f$D)
      delta <- last_inner_loss - inner_loss
      
      # check convergence
      if(delta <= inner_tol){
        break       
      } else{
        last_inner_loss <- inner_loss  
      }
    }
    
    ########## prior variance update ########## 
    for(q in 1:n_factors){
      f <- prior_variance_update(f, q)
    }
    
    ########## residual precision update ########## 
    f <- residual_precision_update(f)
    
    ########## mean update ########## 
    f <- mean_update(f)
    
    ########## check convergence ########## 
    loss[i] <- compute_loss(f$tau, f$S, f$I_n, f$L, f$D)
    delta <- loss[i-1] - loss[i] # should be positve if loss decreasing
    
    # print an update
    if(i %% n_print == 0){
      msg <- paste0("iteration=", i, " | loss=", loss[i], 
                    " | delta=", delta, " | m=", m, " | last_t=", f$t)
      print(msg)
    }
    
    # check convergence
    if(delta <= tol){
      f <- add_convergence_info(f, loss, i)
      return(f)
    }
  }
  
  f <- add_convergence_info(f, loss, i)
  return(f)
  
}