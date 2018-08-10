#' Gotta think of a good description for this...
#'
#' @param theta mean of spike slab
#' @param y value of y
#' @param sigma standard deviation of sampling distribution
#' @param tau standard deviation of slab normal distribution
#' @param t truncation point
#'
int_fun_gauss <- function(theta, y, sigma, tau, t) {
  dnorm(y, theta, sigma) * dnorm(theta, 0, tau) /
    Pr_S_conditional(theta, sigma, t)
}
