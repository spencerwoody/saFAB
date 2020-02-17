#' Density of marginal distribution of y ~ N(theta, sigma^2),
#' where theta has spike-slab prior, under joint selection
#'
#' @param y value at which to evaluate
#' @param sigma standard deviation of sampling distribution
#' @param t truncation point
#' @param p proportion for spike part
#' @param tau standard deviation of normal
#'
#' @export
marginal_ss_joint1 <- function(y, sigma, t, p, mu, tau) {

  density <- (p * dnorm(y, mu, sqrt(sigma^2 + tau^2)) +
                (1 - p) * dnorm(y, 0, sigma) ) *
    ifelse(y > t, 1, 0)

  return(density)

}
