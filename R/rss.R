##' Generate n samples from the spike-slab normal prior, a mixture of
##' a normal slab and point mass spike at zero
##'
##' theta ~ p * 1(theta = 0) + (1 - p) * N(theta; 0, tau^2)
##' @title rss
##' @param n number of samples to generate
##' @param p proportion of zeros
##' @param tau standard deviation of normal (slab) prior
##' @return n samples from the spike-and-slab prior
##'
##' @export
##' @author Spencer Woody
rss <- function(n, p, tau) {
  rnorm(n, 0, tau) * rbinom(n, 1, p)
}
