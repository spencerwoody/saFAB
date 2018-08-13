
bayes_pdf_signal <- function(theta, y, p, sigma, tau) {
  
  v2 <- 1 / (sigma^-2 + tau^-2)
  a <- v2 / sigma^2
  
  dnorm(theta, a * y, sqrt(v2))
  
}

bayes_cdf_quantile <- function(theta, y, p, sigma, tau, q = 0) {
  
  phat_y <- phat(y, p, sigma, tau = tau)
  
  cdf <- phat_y * integrate(bayes_pdf_signal, lower = -Inf, upper = theta,
                            y = y, p = p, sigma = sigma, tau = tau)$value + 
    (1 - phat_y) * (theta >= 0)
  
  return(cdf - q)
  
}

##' Find quantile of Bayesian posterior distribution for normal
##' sampling model, and spike-and-slab prior for theta
##'
##' 
##' @title bayes_find_quantile
##' @param y value of y
##' @param p prior probability of Gaussian slab signal
##' @param sigma standard deviation of sampling distribution
##' @param tau standard deviation / scale of Gaussian slab
##' @param q desired quartile
##' @return quartile from Bayesian posterior distribution
##' @author Spencer Woody
##'
##' @export
bayes_find_quantile <- function(y, p, sigma, tau, q) {
  # q is the desired quantile
  
  root <- uniroot(bayes_cdf_quantile, 
                  lower = min(y - 7, -1), 
                  upper = max(y + 7, 1),
                  y = y, p = p, sigma = sigma, tau = tau, q = q)$root
  
  if (abs(root) < 1e-4) {
    return(0)
  } else {
    return(root)
  }
  
}

