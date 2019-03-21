
## ------------------------------------------------------------------------
## Gaussian prior for theta

##' Marginal density of y for Gaussian prior on mean theta, under
##' joint selection
##'
##' 
##' @title marginal_gauss_joint
##' @param y value of y
##' @param sigma standard deviation
##' @param t truncation point
##' @param tau scale / standard deviation of Gaussian prior
##' @return Marginal density of y for Gaussian prior on mean theta
##' @author Spencer Woody
##'
##' @export
marginal_gauss_joint <- function(y, sigma, t, tau) {
  dnorm(y, 0, sqrt(sigma^2 + tau^2)) * ifelse(abs(y) > t, 1, 0)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param y 
##' @param sigma 
##' @param t 
##' @param mu 
##' @param tau 
##' @return 
##' @author Spencer Woody
##'
##' @export
marginal_gaussnc_joint <- function(y, sigma, t, mu, tau) {
  dnorm(y, mu, sqrt(sigma^2 + tau^2)) * ifelse(abs(y) > t, 1, 0)
}


## ------------------------------------------------------------------------
## Laplace prior for theta

pdf_laplace <- function(theta, tau) {
  # PDF of Laplace distribution, theta ~ Laplace(0, tau) (scale is tau)
  
  exp(- abs(theta) / tau) / (2 * tau)
  
}


int_fun_laplace_joint <- function(theta, y, sigma, tau, t) {
  dnorm(y, theta, sigma) * pdf_laplace(theta, tau)
}

##' marginal density of y for Laplacian prior on mean theta, under
##' joint selection
##'
##' 
##' @title marginal_laplace_joint
##' @param y value of y
##' @param sigma standard deviation for sampling model
##' @param t truncation point
##' @param tau scale for laplacian prior of theta
##' @return density of marginal of y for laplacian prior
##' @author Spencer Woody
##'
##' @export
marginal_laplace_joint <- function(y, sigma, t, tau) {
  integrate(int_fun_laplace_joint,
            lower = -Inf,
            upper = +Inf,
            y = y, sigma = sigma, tau = tau, t = t)$value * 
    ifelse(abs(y) > t, 1, 0)
}



## ------------------------------------------------------------------------
## Horseshoe prior for theta

int_fun_hs_joint <- function(theta, tau, y, sigma) {
  dnorm(y, theta, sigma) * log(1 + 4 * tau^2 / theta^2) 
}

##' marginal density of y for horseshoe prior on mean theta, under
##' joint selection
##'
##' Uses lower bound of horseshoe density, from Carvahlo et al. (2008)
##' @title marginal_hs_joint
##' @param y value of y
##' @param sigma standard deviation of sampling distribution
##' @param tau (global) scale of horseshoe prior
##' @param t truncation point
##' @return marginal density of y for horseshoe prior on mean theta
##' @author Spencer Woody
##'
##' @export
marginal_hs_joint <- function(y, sigma, tau, t) {
  integrate(int_fun_hs_joint, lower = -Inf, upper = Inf, 
            tau = tau, y = y, sigma = sigma)$value * 
    ifelse(abs(y) >= t, 1, 0)
}
