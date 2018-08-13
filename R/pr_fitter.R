
# Function used for estimating marginal of y
int_fun_est <- function(theta, y, sig0, prior_fit, t) {
  dnorm(y, theta, sigma) * prior_fit(theta) / Pr_S_conditional(theta, sig0, t)
}


marginal_est_fixed <- function(y, pi0, sigma, prior_fit, t) {
  
  p <- 1 - pi0
  
  int_part <- integrate(int_fun_est, 
                        lower = -Inf, 
                        upper = +Inf, 
                        y = y, sigma = sigma, prior_fit = prior_fit, t = t)$value
  
  density <- (p * int_part + 
                (1 - p) * dnorm(y, 0, sigma) / Pr_S_fixed(0, sigma, t)) * 
    ifelse(abs(y) > t, 1, 0)
  
  return(density)
  
}

#' Use predictive recursion to obtain a semiparametric estimate of
#' marginal of y,
#'
#' @param y sample
#' @param mu0 foobar
#' @param sig0 foobar
#' @param nulltype foobar
#' @return a list, where my_fit is the fitted marginal or y
#'
#' @export
pr_fitter <- function(y, mu0=NULL, sig0=NULL, nulltype = 'theoretical') {
  require(FDRreg)
  require(splines)

  # sampling model of y
  if(nulltype == 'theoretical') {
    if(missing(mu0)) mu0 = 0
    if(missing(sig0)) sig0 = 1
  } else if(nulltype == 'empirical') {
    null_fit <- efron(y, nmids = 200, df = 15,
                     nulltype = 'empirical')
    mu0 <- null_fit$mu0
    sig0 <- null_fit$sig0
  } else {
    error('nulltype must be either "theoretical" or "empirical"')
  }

  # fit predictive recursion
  pr_fit <- prfdr(y, mu0, sig0)

  # interpolate on the log scale to get the prior
  logprior_fit <- splinefun(pr_fit$x_grid,
                            log(pr_fit$pitheta_grid),
                            method='natural')
  prior_fit <- function(theta) exp(logprior_fit(theta))

  # interpolate on the log scale to get the marginal density
  logmy_fit <- splinefun(pr_fit$x_grid, log(pr_fit$fmix_grid), method='natural')

  ## Marginal for data, with no selection
  
  my_fit <- function(y) exp(logmy_fit(y))

  ## Marginal under joint selection

  my_fit_joint <- function(y, sigma, t) {
    my_fit(y) * ifelse(abs(y) > t, 1, 0)
  }

  ## Marginal under conditional selection

  p <- 1 - pr_fit$pi0

  ## sigma argument needed to match up with Hprime_w_safab function...
  my_fit_conditional <- function(y, sigma, t) {
    int_part <- integrate(int_fun_est, lower = -Inf, upper = +Inf,
                          y = y, sig0 = sig0, prior_fit = prior_fit,
                          t = t)$value
    
    density <- (p * int_part + 
                (1 - p) * dnorm(y, 0, sig0) / Pr_S_conditional(0, sig0, t)) * 
      ifelse(abs(y) > t, 1, 0)
    
    return(density)
  }
  
  return(list(
    prior_fit = prior_fit,
    my_fit = my_fit,
    my_fit_joint = my_fit_joint,
    my_fit_conditional = my_fit_conditional,
    pr_fit = pr_fit,
    sig0 = sig0
  ))

}
