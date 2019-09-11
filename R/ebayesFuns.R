
logit <- function(p) {
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}


marginal_ss <- function(y, sigma, p, tau) {
  ## Marginal of y using spike-slab-prior, but withOUT selection
  
  p * dnorm(y, 0, sqrt(sigma^2 + tau^2)) + 
    (1 - p) * dnorm(y, 0, sigma)
  
}

marginal_y_optim <- function(pars, y) {
  ## Wrapper function for feeding into optim
  
  ## Inputs of par (which are in the set of reals)
  ## pars[1] = log(sigma)
  ## pars[2] = logit(p)
  ## pars[3] = log(tau)

  sigma <- exp(pars[1])
  p <- inv_logit(pars[2])
  tau <- exp(pars[3])
  
  ## Negative log likelihood
  - sum(log(marginal_ss(y, sigma, p, tau)))
  
}

marginal_y_optim_knownsigma <- function(pars, sigma, y) {
  ## Wrapper function for feeding into optim
  
  ## Inputs of par (which are in the set of reals)
  ## pars[1] = log(sigma)
  ## pars[2] = logit(p)
  ## pars[3] = log(tau)

  p <- inv_logit(pars[1])
  tau <- exp(pars[2])
  
  ## Negative log likelihood
  - sum(log(marginal_ss(y, sigma, p, tau)))
  
}

##' Find empirical Bayes estimates of hyperparameters of spike-slab
##' prior by using marginal likelihood maximization
##'
##' Theta follows the spike-slab prior
##' @title eb_fun
##' @param y vector of oberservation
##' @return MMLE estimates of hyperparameters, the hessian matrix at
##'   the MMLE estimate, and the convergence message
##' @author Spencer Woody
##'
##' @export
eb_fun <- function(y) {
  
  eb_optim <- optim(par = c(log(5), log(0.5), log(5)),
                    fn = marginal_y_optim, 
                    y = y,
                    hessian = T)
  
  return(
    list(
      sigma_hat = eb_optim$par[1] %>% exp(),
      p_hat = eb_optim$par[2] %>% inv_logit(),
      tau_hat = eb_optim$par[3] %>% exp(),
      H = eb_optim$hessian,
      covergence = eb_optim$convergence
    )
  )
  
}

##' Find empirical Bayes estimates of hyperparameters of spike-slab
##' prior by using marginal likelihood maximization
##'
##' Theta follows the spike-slab prior
##' @title eb_fun
##' @param y vector of oberservation
##' @return MMLE estimates of hyperparameters, the hessian matrix at
##'   the MMLE estimate, and the convergence message
##' @author Spencer Woody
##'
##' @export
eb_fun_knownsigma <- function(y, sigma) {
  
  eb_optim <- optim(par = c(log(0.5), log(5)),
                    fn = marginal_y_optim_knownsigma, 
                    y = y, sigma = sigma,
                    hessian = T)
  
  return(
    list(
      p_hat = eb_optim$par[1] %>% inv_logit(),
      tau_hat = eb_optim$par[2] %>% exp(),
      H = eb_optim$hessian,
      covergence = eb_optim$convergence
    )
  )
  
}
