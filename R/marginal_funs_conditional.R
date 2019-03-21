
## ------------------------------------------------------------------------
## Gaussian prior for theta

##' Marginal density of y for Gaussian prior on mean theta, under
##' conditional selection
##'
##' 
##' @title marginal_gauss_conditional
##' @param y value of y
##' @param sigma standard deviation of sampling distribution
##' @param tau scale / standard deviation of Gaussian prior for theta
##' @param t truncation point
##' @return marginal density
##' @author Spencer Woody
##'
##' @export
marginal_gauss_conditional <- function(y, sigma, tau, t) {
  integrate(int_fun_gauss,
            lower = -Inf, 
            upper = +Inf,
            y = y, sigma = sigma, tau = tau, t = t)$value * 
    ifelse(abs(y) > t, 1, 0)
}


marginal_gaussnc_conditional <- function(y, sigma, mu, tau, t) {
  integrate(int_fun_gaussnc,
            lower = -Inf, 
            upper = +Inf,
            y = y, sigma = sigma, mu = mu, tau = tau, t = t)$value * 
    ifelse(abs(y) > t, 1, 0)
}



## ------------------------------------------------------------------------
## Laplace prior for theta

int_fun_laplace_conditional <- function(theta, y, sigma, tau, t) {
  dnorm(y, theta, sigma) * pdf_laplace(theta, tau) / Pr_S_conditional(theta, sigma, t)
}



##' Marginal density of y for Laplacian prior on mean theta, under
##' conditional selection
##'
##' 
##' @title marginal_laplace_conditional
##' @param y value of y
##' @param sigma standard deviation of sampling distribution
##' @param tau scale of Laplacian prior for theta
##' @param t truncation point
##' @return marginal density
##' @author Spencer Woody
##'
##' @export
marginal_laplace_conditional <- function(y, sigma, tau, t) {
  integrate(int_fun_laplace_conditional,
            lower = -Inf,
            upper = +Inf,
            y = y, sigma = sigma, tau = tau, t = t)$value * 
    ifelse(abs(y) > t, 1, 0)
}

## ------------------------------------------------------------------------
## Horseshoe prior for theta

int_fun_hs_conditional <- function(theta, tau, y, sigma, t) {
  dnorm(y, theta, sigma) * log(1 + 4 * tau^2 / theta^2) / 
    Pr_S_conditional(theta, sigma, t)
}


##' Marginal density of y for horseshoe prior on mean theta, under
##' conditional selection
##'
##' Uses closed-form lower bound density of horseshoe, as described by
##' Carvalho et al. (2008)
##' @title marginal_hs_conditional
##' @param y value of y
##' @param sigma standard deviation of sampling distribution
##' @param tau (global) scale of horse
##' @param t truncation point
##' @return marginal density
##' @author Spencer Woody
##'
##' @export
marginal_hs_conditional <- function(y, sigma, tau, t) {
  integrate(int_fun_hs_conditional, lower = -Inf, upper = Inf, 
            tau = tau, y = y, sigma = sigma, t = t)$value * 
    ifelse(abs(y) >= t, 1, 0)
}

