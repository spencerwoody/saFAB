

## ------------------------------------------------------------------------
## Functions for selection-adjusted Bayesian inference (i.e. under
## conditional selection)

sabayes_pdf_signal <- function(theta, y, p, sigma, tau, t) {
  # Selection-adjusted posterior of theta, assuming it's signal
  # \pi_S(\theta | y, \gamma = 1)
  
  v2 <- 1 / (sigma^-2 + tau^-2)
  a <- v2 / sigma^2
  
  dnorm(theta, a * y, sqrt(v2)) / 
    (pnorm((-t - theta) / sigma) + 1 - pnorm((t - theta) / sigma))
}


make_d <- function(y, p, sigma, tau, t) {
  # Normalizing constant for saposterior_pdf_signal
  
  integrate(sabayes_pdf_signal, lower = -Inf, upper = +Inf, 
            y = y, p = p, sigma = sigma, tau = tau, t = t)$value
  
}

sabayes_cdf_quantile <- function(theta, y, p, sigma, tau, t, d = NULL,
                                 q = 0) {
  # Calculates selection-adjusted posterior CDF, with possible subtraction 
  # of desired quantile for root-finding
  
  phat_y <- phat(y, p, sigma, tau = tau)
  
  if (is.null(d)) {
    d <- make_d(y, p, sigma, tau, t)
  }
  
  weight_spike <- (1 - phat_y) / (pnorm(-t / sigma) + 1 - pnorm(t / sigma))
  weight_slab <- d * phat_y
  
  myint <- integrate(sabayes_pdf_signal, lower = -Inf, upper = theta,
                     y = y, p = p, sigma = sigma, tau = tau, t = t)$value
  
  cdf <- weight_spike / (weight_spike + weight_slab) * (theta >= 0) + 
    (weight_slab / (weight_spike + weight_slab) * myint / d)
  
  return(cdf - q)
  
}

##' Find quantile of selection-adjusted Bayesian posterior
##' distribution (i.e. under conditional selection), using normal
##' sampling distribution and spike-and-slab prior on means thetas
##'
##' .. content for \details{} ..
##' @title sabayes_find_quantile
##' @param y value of y
##' @param p prior probability of signal
##' @param sigma standard deviation of sampling distribution
##' @param tau standard deviation / scale of Gaussian slab
##' @param t truncation point
##' @param q desired quantile
##' @return Quantile of Bayesian posterior distribution
##' @author Spencer Woody
##'
##' @export
sabayes_find_quantile <- function(y, p, sigma, tau, t, q) {
  # q is the desired quantile
  
  d <- make_d(y, p, sigma, tau, t)
  
  root <- uniroot(sabayes_cdf_quantile, 
                  lower = min(y - 5, -5), upper = max(y + 5, 5),
                  y = y, p = p, sigma = sigma, tau = tau, t = t, d = d, q = q)$root
  
  if (abs(root) < 1e-4) {
    return(0)
  } else {
    return(root)
  }
  
}

