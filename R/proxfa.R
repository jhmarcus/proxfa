

#' @title Proximal gradient descent for factor analysis
#'
#' @description This function computes maximum likelihood estimates of the parameters in a 
#'  factor analysis model with constraints on the loadings. The model can be motivated 
#'  by the following generative factor analysis model:
#'  
#'    f_j ~ N(0, D)
#'    y_j | f_j ~ N(\mu_j1 + Lf_j, \tau^{-1}I)
#'  
#'  Here y_j is the data, f_j are called latent factors and L are sample loadings on to these factors. 
#'  We consider the case where L has some constraint, such as non-negativity, and our goal is to 
#'  estimate L by maximum likelihood, in addition to other fixed parameters. Typically models 
#'  like this are fit fit using an EM algorithim. Here we take a direct approach and marginalize 
#'  out the latent factors:
#'  
#'    y_j ~ N(\mu_j1, LDL^T + \tau^{-1}I) 
#'  
#'  We then optimize marginal likelihood using alternating optimization and proximal gradient descent.
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
#' @param fix_step_size is a boolean to fix the step size rather than line search
#' 
#' @return returns a list with fitted values
#' 
#' @export
#' 
proxfa <- function(Y, L, D, mu, tau, 
                   step_size, step_size_shrink=.5,
                   max_iter=500, max_inner_iter=20, 
                   tol=1e-3, inner_tol=1e-3, n_print=10, 
                   prox_fn=identity_prox_fn, eps=1e-8, fix_D=FALSE, 
                   fix_mu=FALSE, fix_tau=FALSE, fix_step_size=FALSE){
  
  n <- nrow(Y) # number of samples
  p <- ncol(Y) # number of features
  n_factors <- ncol(L) # number of factors
  
  if(n_factors <=2 && fix_D ==FALSE){
    stop("D updates not implemented for Q <= 2")
  }
  
  ones_n <- rep(1, n)
  
  # compute sample covarinace including the feature means
  YYt <- Y %*% t(Y)
  S <- comp_samp_cov_wmean(YYt, Y, mu, ones_n)
  
  # intial fit object
  f <- list(Y=Y, YYt=YYt, n=n, p=p, 
            step_size=step_size, step_size_shrink=step_size_shrink, 
            fix_step_size=fix_step_size, t=NA,
            S=S, L=L, D=D, mu=mu, tau=tau, 
            ones_n=ones_n, I_n=diag(ones_n), 
            prox_fn=prox_fn, eps=eps)
      
  # setup loss      
  loss = rep(NA, max_iter+1)
  loss[1] <- comp_neg_loglik(f$tau, f$S, f$I_n, f$L, f$D)
  for(i in 2:(max_iter+1)){
    
    ########## loadings update ########## 
    last_inner_loss <- loss[1]
    
    # gradient descent
    for(m in 1:max_inner_iter){
      f <- loadings_update(f)  
      inner_loss <- comp_neg_loglik(f$tau, f$S, f$I_n, f$L, f$D)
      delta <- last_inner_loss - inner_loss
      
      # check convergence
      if(delta <= inner_tol){
        break       
      } else{
        last_inner_loss <- inner_loss  
      }
    }
    
    ########## prior variance update ########## 
    if(!fix_D){
      for(q in 1:n_factors){
        f <- prior_variance_update(f, q)
      }      
    }

    ########## residual precision update ########## 
    if(!fix_tau){
      f <- residual_precision_update(f)
    }

    ########## mean update ########## 
    if(!fix_mu){
      f <- mean_update(f)
    }
    
    ########## check convergence ########## 
    loss[i] <- comp_neg_loglik(f$tau, f$S, f$I_n, f$L, f$D)
    delta <- loss[i-1] - loss[i] # should be positve if loss decreasing
    
    # print an update
    if(i %% n_print == 0){
      msg <- paste0("iteration=", i, 
                    " | loss=", loss[i], 
                    " | delta=", delta, 
                    " | m=", m, 
                    " | last_t=", f$t)
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