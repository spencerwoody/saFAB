
## Non-selective spending functions
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param y 
##' @param sigma 
##' @param mu 
##' @param tau 
##' @return 
##' @author Spencer Woody
##'
##' @export
marginal_gauss <- function(y, sigma, mu, tau) {
  dnorm(y, mu, sqrt(sigma^2 + tau^2)) 
}

## marginalCdf_gauss <- function(y, sigma, )
