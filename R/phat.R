##' Posterior probability of being a signal, using spike-slab prior
##'
##' Mixture of point-mass spike at zero, and Gaussian slab
##' @title phat
##' @param y value of y
##' @param p prior probability of signal (slab)
##' @param sigma standard deviation of sampling distribution of y
##' @param tau standard deviation / scale of Gaussian slab
##' @return posterior probability of signal
##' @author Spencer Woody
phat <- function(y, p, sigma, tau) {
  p * dnorm(y, 0, sqrt(sigma^2 + tau^2)) / 
    (p * dnorm(y, 0, sqrt(sigma^2 + tau^2)) + (1 - p) * dnorm(y, 0, sigma))
}
