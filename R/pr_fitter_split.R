
# Function used for estimating marginal of y
int_fun_cnd <- function(theta, y, sigma_orig, n = 1, prior_fit, t) {
  dnorm(y, theta, sigma_orig / sqrt(n)) * prior_fit(theta) /
    Pr_S_cnd(theta, sigma_orig / sqrt(n), t)
}

int_fun_jnt <- function(theta, y, sigma_orig, n = 1, prior_fit, t) {
  dnorm(y, theta, sigma_orig / sqrt(n)) * prior_fit(theta) 
}

##' Use predictive recursion to obtain a semiparametric estimate of
##' marginal of y,
##'
##' @param y sample
##' @param mu0 foobar
##' @param sig0 foobar
##' @param nulltype foobar
##' @return a list, where my_fit is the fitted marginal or y
##'
##' @export
pr_fitter_splitmat <- function(y, mu0=NULL, sig0=NULL, nulltype = "theoretical") {
  require(FDRreg)
  require(splines)

  ## Split y's
  N <- nrow(y)
  yTr <- y[, 1] #y's for estimating prior

  ## y's for constructing intervals
  yInt <- matrix(y[, -1], nrow = N) #y's for constructing intervals

  n <- ncol(yInt)

  ## sampling model of y
  if(nulltype == "theoretical") {
    if(missing(mu0)) mu0 = 0
    if(missing(sig0)) sig0 = 1
  } else if(nulltype == "empirical") {
    null_fit <- efron(yTr, nmids = 200, df = 15,
                     nulltype = "empirical")
    mu0 <- null_fit$mu0
    sig0 <- null_fit$sig0
  } else {
    error("nulltype must be either \"theoretical\" or \"empirical\"")
  }

  ## fit predictive recursion
  pr_fit <- prfdr(yTr, mu0, sig0)

  ## interpolate on the log scale to get the prior
  logprior_fit <- splinefun(pr_fit$x_grid,
                            log(pr_fit$pitheta_grid),
                            method='natural')
  prior_fit <- function(theta) exp(logprior_fit(theta))

  ## interpolate on the log scale to get the marginal density
  logmy_fit <- splinefun(pr_fit$x_grid, log(pr_fit$fmix_grid), method="natural")

  ## Marginal for data, with no selection
  
  my_fit <- function(y) exp(logmy_fit(y))

  ## Marginal under joint selection

  p <- 1 - pr_fit$pi0

  my_fit_jnt <- function(y, sigma, t, sigma_orig, n = n) {

    int_part <- integrate(int_fun_jnt, 
                          lower = -Inf, 
                          upper = +Inf, n = n,
                          y = y, sigma_orig = sigma_orig,
                          prior_fit = prior_fit, t = t)$value

    (p * int_part + (1 - p) * dnorm(y, 0, sigma_orig / sqrt(n))) *
    ifelse(abs(y) > t, 1, 0)
    
  }
  
  ## Marginal under conditional selection



  ## sigma argument needed to match up with Hprime_w_safab function...
  my_fit_cnd <- function(y, sigma, t, sigma_orig, n = n) {
    int_part <- integrate(int_fun_cnd, lower = -Inf, upper = +Inf,
                          y = y, sigma_orig = sigma_orig,
                          prior_fit = prior_fit, n = n, t = t)$value
    
    (p * int_part + 
     (1 - p) * dnorm(y, 0, sigma_orig / sqrt(n)) /
     Pr_S_cnd(0, sigma_orig / sqrt(n), t)) * 
      ifelse(abs(y) > t, 1, 0)
    
  }
  
  return(list(
    prior_fit = prior_fit,
    my_fit = my_fit,
    my_fit_jnt = my_fit_jnt,
    my_fit_cnd = my_fit_cnd,
    pr_fit = pr_fit,
    sig0 = sig0,
    phat = p
  ))

}




marginal_np_cnd <- function(y, pi0, sigma, n = 1, prior_fit, t) {
  
  p <- 1 - pi0
  
  int_part <- integrate(int_fun_cnd, 
                        lower = -Inf, 
                        upper = +Inf, 
                        y = y, sigma = sigma, n = n,
                        prior_fit = prior_fit, t = t)$value
  
  (p * int_part + 
   (1 - p) * dnorm(y, 0, sigma / sqrt(n)) /
   Pr_S_cnd(0, sigma / sqrt(n), t)) * ifelse(abs(y) > t, 1, 0)
  
  
}

marginal_np_jnt <- function(y, pi0, sigma, n = 1, prior_fit, t) {

  p <- 1 - pi0
  
  int_part <- integrate(int_fun_jnt, 
                        lower = -Inf, 
                        upper = +Inf, 
                        y = y, sigma = sigma, prior_fit = prior_fit, t = t)$value

  (p * int_part + (1 - p) * dnorm(y, 0, sigma / sqrt(n))) * ifelse(abs(y) > t, 1, 0)
  
}


