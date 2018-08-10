#' Density of marginal distribution of y ~ N(theta, sigma^2),
#' where theta has spike-slab prior, under joint selection
#'
#' @param y value at which to evaluate
#' @param sigma standard deviation of sampling distribution
#' @param t truncation point
#' @param p proportion for spike part
#' @param tau standard deviation of normal
#'
marginal_ss_joint <- function(y, sigma, t, p, tau) {

  density <- (p * dnorm(y, 0, sqrt(sigma^2 + tau^2)) +
                (1 - p) * dnorm(y, 0, sigma) ) *
    ifelse(abs(y) > t, 1, 0)

  return(density)

}
